clear all;
close all;
clc;

%% Przetwarzanie wszystkich obrazów w określonym folderze
folder_path = 'D:\Biometria\Iris2\Integrodifferential operator\CSIL';
output_folder = 'D:\Biometria\Iris2\Integrodifferential operator\database';

% Pobierz listę plików obrazów z katalogu OFTA
image_files = dir(fullfile(folder_path, '*.jpg'));

% Pętla po każdym pliku obrazu
for k = 1:length(image_files)
    % Wczytaj obraz
    img_filename = fullfile(folder_path, image_files(k).name);
    img_org = imread(img_filename);
    img_resized=imresize(img_org,1/2);
    
    % Przetwarzanie wstępne
    img_gray = im2gray(img_resized); % Konwersja obrazu na skalę szarości
    img = double(img_gray); % Konwersja na typ double

    % Detekcja krawędzi za pomocą operatora Sobela
    sx = [-1 0 1; -2 0 2; -1 0 1];
    sy = [-1 -2 -1; 0 0 0; 1 2 1];
    grad = zeros(size(img));
    for i = 1:size(img, 1) - 2
        for j = 1:size(img, 2) - 2
            Gx = sum(sum(sx .* img(i:i+2, j:j+2))); % Gradient w kierunku X
            Gy = sum(sum(sy .* img(i:i+2, j:j+2))); % Gradient w kierunku Y
            grad(i+1, j+1) = sqrt(Gx.^2 + Gy.^2); % Łączny gradient
        end
    end
    threshold = 30;
    grad = uint8(grad);
    grad = grad > threshold; % Próg detekcji krawędzi

    %% Transformata Hougha
    [rows, cols] = size(grad);

    % Inicjalizacja parametrów
    min_iris_radius = 31;
    max_iris_radius = 90;
    min_pupil_radius = 15;
    max_pupil_radius = 30;
    iris_threshold = 0.7;
    pupil_threshold = 0.7;

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

    %% Wycinanie tęczówki
    theta = 0:0.01:2*pi;
    radii = linspace(iris_r, pupil_r, 128); % Generowanie promieni
    [thetaGrid, radiiGrid] = meshgrid(theta, radii);
    % Konwersja współrzędnych biegunowych na kartezjańskie
    x = iris_x + radiiGrid .* cos(thetaGrid);
    y = iris_y + radiiGrid .* sin(thetaGrid);
    % Normalizacja obrazu tęczówki za pomocą griddata
    normalizedIris = interp2(double(img), x, y, 'linear', 0);
    % Normalizacja rozwiniętej tęczówki do stałych wymiarów
    
   

   % Parametry filtra Gabora
wavelength = 8;      % Długość fali
theta = 0;    % Orientacja w radianach
sigmaX = 4;          % Standardowe odchylenie w osi X
sigmaY = 4;          % Standardowe odchylenie w osi Y
phaseOffset = 0;     % Przesunięcie fazowe

% Definicja siatki współrzędnych
[x, y] = meshgrid(-fix(sigmaX*3):fix(sigmaX*3), -fix(sigmaY*3):fix(sigmaY*3)); %generacja

    
% Obrót współrzędnych
xPrime = x .* cos(theta) + y .* sin(theta);
yPrime = -x .* sin(theta) + y .* cos(theta);
    
% Generowanie rzeczywistej i urojonej części filtra
realPart = exp(-0.5 * (xPrime.^2 / sigmaX^2 + yPrime.^2 / sigmaY^2)) .* ...
               cos(2 * pi * xPrime / wavelength + phaseOffset);
imagPart = exp(-0.5 * (xPrime.^2 / sigmaX^2 + yPrime.^2 / sigmaY^2)) .* ...
               sin(2 * pi * xPrime / wavelength + phaseOffset);
    
% Tworzenie filtra zespolonego
gaborFilter = realPart + 1i * imagPart;



% Get the size of the image
[height, width] = size(normalizedIris);

% Calculate the height and width of each segment
segmentHeight = height / 16;
segmentWidth = (width-1) / 157 ;

% Initialize a 3D array to store the segments
segments = zeros(segmentHeight, segmentWidth, 16*157);

% Initialize an array to store the mean values
meanValues = zeros(1, 16*157);

% Loop over each segment
for i = 1:16
    for j = 1:157
        % Calculate the start and end rows and columns for this segment
        startRow = round((i - 1) * segmentHeight) + 1;
        endRow = round(i * segmentHeight);
        startCol = round((j - 1) * segmentWidth) + 1;
        endCol = round(j * segmentWidth);

        % Extract the segment from the image
        segment = normalizedIris(startRow:endRow, startCol:endCol);

        % Store the segment
        segments(:,:,((i-1)*157)+j) = segment;

        % Calculate the mean of the segment and store it
        meanValues(((i-1)*157)+j) = mean(segment(:));
    end
end



realFilteredImg = filter2(real(gaborFilter), meanValues);
imagFilteredImg = filter2(imag(gaborFilter), meanValues);
filteredImgComplex = realFilteredImg + 1i * imagFilteredImg;      
% Assuming 'filteredImgComplex' is your filtered image
% Binarize the image
threshold = mean(filteredImgComplex(:));
binary_image = filteredImgComplex > threshold;
figure
% Encode the binary image into a bit string
iris_code = reshape(binary_image, 1, []);

    % Zapisz wynik do pliku .mat
    [~, name, ~] = fileparts(img_filename);
    output_filename = fullfile(output_folder, [name, '_iris_code.mat']);
    iris_code_variable_name = [name, '_iris_code'];
    eval([iris_code_variable_name, ' = iris_code;']);
    save(output_filename, iris_code_variable_name);

    disp(['Saved iris code to ', output_filename]);
end
