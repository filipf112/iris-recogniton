clear all;
close all;
clc;

%% Przetwarzanie wszystkich obrazów w określonym folderze
folder_path = 'D:\Biometria\Iris2\Integrodifferential operator\OFTA';
output_folder = 'D:\Biometria\Iris2\Integrodifferential operator\database1';

% Pobierz listę plików obrazów z katalogu OFTA
image_files = dir(fullfile(folder_path, '*.bmp'));

% Pętla po każdym pliku obrazu
for k = 1:length(image_files)
    % Wczytaj obraz
    img_filename = fullfile(folder_path, image_files(k).name);
    img_org = imread(img_filename);
    img_height = 300;
    img_width = 444;
    img_resized = imresize(img_org, [img_height, img_width]);
    
    % Przetwarzanie wstępne
    img_gray = rgb2gray(img_resized);
    img = double(img_gray);

    % Detekcja krawędzi za pomocą operatora Sobela
    sx = [-1 0 1; -2 0 2; -1 0 1];
    sy = [-1 -2 -1; 0 0 0; 1 2 1];
    grad = zeros(size(img));
    for i = 1:size(img, 1) - 2
        for j = 1:size(img, 2) - 2
            Gx = sum(sum(sx .* img(i:i+2, j:j+2)));
            Gy = sum(sum(sy .* img(i:i+2, j:j+2)));
            grad(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
        end
    end
    threshold = 80;
    grad = uint8(grad);
    grad = grad > threshold;

    %% Transformata Hougha
    [rows, cols] = size(grad);

    % Inicjalizacja parametrów
    min_iris_radius = 130;
    max_iris_radius = 170;
    min_pupil_radius = 30;
    max_pupil_radius = 80;

    % Inicjalizacja akumulatorów Hougha dla tęczówki i źrenicy
    iris_accumulator = zeros(rows, cols, max_iris_radius);
    pupil_accumulator = zeros(rows, cols, max_pupil_radius);

    % Detekcja tęczówki
    for y = 1:rows
        for x = 1:cols
            if grad(y, x) == 1
                for r = min_iris_radius:max_iris_radius
                    for theta = 0:pi/100:2*pi
                        a = round(x - r * cos(theta));
                        b = round(y - r * sin(theta));
                        if a > 0 && a <= cols && b > 0 && b <= rows
                            iris_accumulator(b, a, r) = iris_accumulator(b, a, r) + 1;
                        end
                    end
                end
            end
        end
    end

    % Detekcja źrenicy (w ograniczonym obszarze)
    for y = 1:rows
        for x = 1:cols
            if grad(y, x) == 1
                for r = min_pupil_radius:max_pupil_radius
                    for theta = 0:pi/100:2*pi
                        a = round(x - r * cos(theta));
                        b = round(y - r * sin(theta));
                        if a > 0 && a <= cols && b > 0 && b <= rows
                            pupil_accumulator(b, a, r) = pupil_accumulator(b, a, r) + 1;
                        end
                    end
                end
            end
        end
    end
    

    % Znajdź najjaśniejszy punkt w akumulatorze dla tęczówki
    [~, max_iris_index] = max(iris_accumulator(:));
    [iris_y, iris_x, iris_r] = ind2sub(size(iris_accumulator), max_iris_index);

    % Znajdź najjaśniejszy punkt w akumulatorze dla źrenicy
    [~, max_pupil_index] = max(pupil_accumulator(:));
    [pupil_y, pupil_x, pupil_r] = ind2sub(size(pupil_accumulator), max_pupil_index);
    
    figure;
    imshow(img_resized);
    hold on;
    iris_theta = 0:0.01:(2*pi);
    pupil_theta = 0:0.01:(2*pi);
    plot(iris_x + iris_r * cos(iris_theta), iris_y + iris_r * sin(iris_theta), 'b', 'LineWidth', 2);
    plot(pupil_x + pupil_r * cos(pupil_theta), pupil_y + pupil_r * sin(pupil_theta), 'b', 'LineWidth', 2);
    title('Wykryte okręgi tęczówki i źrenicy');
    hold off;
    %% Wycinanie tęczówki
    theta = 0:0.01:2*pi;
    radii = linspace(iris_r, pupil_r, iris_r + pupil_r);
    [thetaGrid, radiiGrid] = meshgrid(theta, radii);
    x = iris_x + radiiGrid .* cos(thetaGrid);
    y = iris_y + radiiGrid .* sin(thetaGrid);
    normalizedIris = interp2(double(img), x, y, 'linear', 0);

    % Normalizacja rozwiniętej tęczówki do stałych wymiarów
    normalized_height = 160;
    normalized_width = 624;
    normalizedIris = imresize(normalizedIris, [normalized_height, normalized_width]);

    %% Gabor Filter
    wavelength = 8;
    theta = 0;
    sigmaX = 4;
    sigmaY = 4;
    phaseOffset = 0;

    [x, y] = meshgrid(-fix(sigmaX*3):fix(sigmaX*3), -fix(sigmaY*3):fix(sigmaY*3));
    xPrime = x .* cos(theta) + y .* sin(theta);
    yPrime = -x .* sin(theta) + y .* cos(theta);
    realPart = exp(-0.5 * (xPrime.^2 / sigmaX^2 + yPrime.^2 / sigmaY^2)) .* cos(2 * pi * xPrime / wavelength + phaseOffset);
    imagPart = exp(-0.5 * (xPrime.^2 / sigmaX^2 + yPrime.^2 / sigmaY^2)) .* sin(2 * pi * xPrime / wavelength + phaseOffset);
    gaborFilter = realPart + 1i * imagPart;     

    [height, width] = size(normalizedIris);
    segmentHeight = height / 16;
    segmentWidth = width / 156;

    segments = zeros(segmentHeight, segmentWidth, 16*156);
    meanValues = zeros(1, 16*156);

    for i = 1:16
        for j = 1:156
            startRow = round((i - 1) * segmentHeight) + 1;
            endRow = round(i * segmentHeight);
            startCol = round((j - 1) * segmentWidth) + 1;
            endCol = round(j * segmentWidth);
            segment = normalizedIris(startRow:endRow, startCol:endCol);
            segments(:,:,((i-1)*156)+j) = segment;
            meanValues(((i-1)*156)+j) = mean(segment(:));
        end
    end

    realFilteredImg = filter2(real(gaborFilter), meanValues);
    imagFilteredImg = filter2(imag(gaborFilter), meanValues);
    filteredImgComplex = realFilteredImg + 1i * imagFilteredImg;      
    
    % Binarize the image
    threshold = mean(filteredImgComplex(:));
    binary_image = filteredImgComplex > threshold;
    iris_code = reshape(binary_image, 1, []);

    % Zapisz wynik do pliku .mat
    [~, name, ~] = fileparts(img_filename);
    output_filename = fullfile(output_folder, [name, '_iris_code.mat']);
    iris_code_variable_name = [name, '_iris_code'];
    eval([iris_code_variable_name, ' = iris_code;']);
    save(output_filename, iris_code_variable_name);

    disp(['Saved iris code to ', output_filename]);
end
