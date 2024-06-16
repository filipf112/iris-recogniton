clear all;
close all;

% Load the fixed iris code (code1)
code1_data = load('D:\Biometria\Iris2\Integrodifferential operator\database1\10o_sr101_iris_code.mat');
varName1 = fieldnames(code1_data);
code1 = code1_data.(varName1{1});

% Directory containing the database of iris codes
database_path = 'D:\Biometria\Iris2\Integrodifferential operator\database1\';
files = dir(fullfile(database_path, '*_iris_code.mat'));

% Path to the OFTA folder
ofta_folder_path = 'D:\Biometria\Iris2\Integrodifferential operator\OFTA';

% Get list of image files in OFTA folder
image_files = dir(fullfile(ofta_folder_path, '*.bmp'));

% Calculate number of files and categorize them
num_files = length(image_files);
authentics = 9;
impostors = num_files - authentics;


% Initialize variables to store FRR and FAR
thresholds = 0.0:0.005:0.5;
FRR = zeros(size(thresholds));
FAR = zeros(size(thresholds));

% Loop through each threshold
for t = 1:length(thresholds)
    threshold = thresholds(t);
    recognized_files = {};
    unrecognized_files = {};
    correct_recognitions = 0;
    false_accepts = 0;
    
    % Loop through each iris code file in the database
    for i = 1:length(files)
        % Load the current iris code (code2)
        code2_data = load(fullfile(database_path, files(i).name));
        varName2 = fieldnames(code2_data);
        code2 = code2_data.(varName2{1});
        
        % Calculate the adaptive Hamming distance with the original code
        distance_original = adaptiveHammingDistance(code1, code2);
        
        % Flip the order of bits in code2
        code2_flipped = [code2(end-511:end), code2(1:end-512)];
        
        % Calculate the adaptive Hamming distance with the flipped code
        distance_flipped = adaptiveHammingDistance(code1, code2_flipped);
        
        % Choose the minimum distance
        distance = min(distance_original, distance_flipped);
        
        % Check if the distance is below the threshold
        if distance <= threshold
            recognized_files{end+1} = files(i).name;
        else
            unrecognized_files{end+1} = files(i).name;
        end
    end
    
    % Extract prefixes from the recognized files
    recognized_prefixes = cellfun(@(x) x(1:6), recognized_files, 'UniformOutput', false);
    
    % Check the number of correctly recognized prefixes
    correct_recognitions = sum(strcmp(recognized_prefixes, '10o_sr'));
    false_accepts = length(recognized_files) - correct_recognitions;
    
end

% Calculate distances between genuine and impostor samples
authentic_distances = []; % distances between genuine samples
impostor_distances = []; % distances between genuine and impostor samples

for i = 1:length(files)
    % Load the current iris code (code2)
    code2_data = load(fullfile(database_path, files(i).name));
    varName2 = fieldnames(code2_data);
    code2 = code2_data.(varName2{1});

    % Calculate the adaptive Hamming distance with the original code
    distance_original = adaptiveHammingDistance(code1, code2);
    
    % Flip the order of bits in code2
    code2_flipped = [code2(end-511:end), code2(1:end-512)];
    
    % Calculate the adaptive Hamming distance with the flipped code
    distance_flipped = adaptiveHammingDistance(code1, code2_flipped);
    
    % Choose the minimum distance
    distance = min(distance_original, distance_flipped);
    
    % Categorize distances based on file type
    if contains(files(i).name, '10o_sr') % genuine sample
        authentic_distances = [authentic_distances, distance];
    else % impostor sample
        impostor_distances = [impostor_distances, distance];
    end
end

% Calculate mean and variance for both authentic and impostor distances
mean_authentic = mean(authentic_distances);
var_authentic = var(authentic_distances);
mean_impostor = mean(impostor_distances);
var_impostor = var(impostor_distances);

% Display means and variances
disp(['Mean Hough for OFTA (authentic): ', num2str(mean_authentic)]);
disp(['Variance Hough for OFTA (authentic): ', num2str(var_authentic)]);
disp(['Mean Hough for OFTA (impostor): ', num2str(mean_impostor)]);
disp(['Variance Hough for OFTA(impostor): ', num2str(var_impostor)]);

% Generate Gaussian PDFs
x = linspace(0, 1, 1000);
pdf_authentic = normpdf(x, mean_authentic, sqrt(var_authentic));
pdf_impostor = normpdf(x, mean_impostor, sqrt(var_impostor));

% Calculate the area under the curves
area_p_a = trapz(x, pdf_authentic);
area_p_i = trapz(x, pdf_impostor);

% Scale the densities
sum_p_a = sum(pdf_impostor);
sum_p_i = sum(pdf_impostor);

p_a_scaled = (pdf_authentic / sum_p_a) ;
p_i_scaled = (pdf_impostor / sum_p_i) ;



% Normalize the densities so that the area under each curve is 1
p_a_normalized = p_a_scaled * area_p_a;
p_i_normalized = p_i_scaled * area_p_i;
% Plot Gaussian PDFs
figure;
hold on;
plot(x, p_a_normalized, 'b', 'LineWidth', 2);
plot(x, p_i_normalized, 'r', 'LineWidth', 2);


xlabel('Odległość Hamminga');
ylabel('Prawdopodobieństwo');
legend('P_{a}(x) (autentics)', 'P_{i}(x) (impostors)');

grid on;
hold off;

% Hamming distance vector
hamming_distance_2 = 0.4:0.005:0.7;

% Initialize FRR and FAR
FRR = zeros(size(hamming_distance_2));
FAR = zeros(size(hamming_distance_2));

% Calculate FRR and FAR for each value in hamming_distance_2
for i = 1:length(hamming_distance_2)
        threshold = hamming_distance_2(i);
    
    % FRR: Area under P_a(x) to the left of the threshold
    FRR(i) = trapz(x(x <= threshold), p_a_normalized(x <= threshold));
    
    % FAR: Area under P_i(x) to the right of the threshold
    FAR(i) = trapz(x(x > threshold), p_i_normalized(x > threshold));
end

% Plot DET curve with logarithmic scale
figure;
hold on;
plot(FAR, FRR, 'b-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('FAR');
ylabel('FRR');
grid on;
hold off;


% Function to calculate the adaptive Hamming distance
function distance = adaptiveHammingDistance(code1, code2)
    % Ensure the codes are the same length
    assert(length(code1) == length(code2), 'Iris codes must be of the same length');
    
    % Calculate the Hamming distance
    distance = sum(code1 ~= code2) / length(code1);
end

