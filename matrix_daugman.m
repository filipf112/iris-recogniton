clear all;
close all;

% Load the fixed iris code (code1)
code1_data = load('D:\Biometria\Iris2\Integrodifferential operator\database1\10o_sr101_iris_code.mat');
varName1 = fieldnames(code1_data);
code1 = code1_data.(varName1{1});

% Directory containing the database of iris codes
database_path = 'D:\Biometria\Iris2\Integrodifferential operator\database1\';
files = dir(fullfile(database_path, '*_iris_code.mat'));

% Path to the CSIL folder
csil_folder_path = 'D:\Biometria\Iris2\Integrodifferential operator\CSIL';

% Get list of image files in CSIL folder
image_files = dir(fullfile(csil_folder_path, '*.jpg'));

% Calculate number of files and categorize them
num_files = length(image_files);
authentics = 9;
impostors = num_files - authentics;

% Initialize variables to store results
threshold = 0.4327;
correct_acceptance = 0;
false_rejection = 0;
false_acceptance = 0;
correct_rejection = 0;

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
    if contains(files(i).name, '10o_sr') % genuine sample
        if distance <= threshold
            correct_acceptance = correct_acceptance + 1;
        else
            false_rejection = false_rejection + 1;
        end
    else % impostor sample
        if distance <= threshold
            false_acceptance = false_acceptance + 1;
        else
            correct_rejection = correct_rejection + 1;
        end
    end
end

% Display the confusion matrix
confusion_matrix = [
    correct_acceptance, false_acceptance;
    false_rejection, correct_rejection
];

figure;
heatmap({'Predicted Authentic', 'Predicted Impostor'}, {'Actual Authentic', 'Actual Impostor'}, confusion_matrix, ...
    'ColorbarVisible', 'on', 'Colormap', parula);
title('Confusion Matrix');
xlabel('Predicted Class');
ylabel('Actual Class');

% Function to calculate the adaptive Hamming distance
function distance = adaptiveHammingDistance(code1, code2)
    % Ensure the codes are the same length
    assert(length(code1) == length(code2), 'Iris codes must be of the same length');
    
    % Calculate the Hamming distance
    distance = sum(code1 ~= code2) / length(code1);
end
