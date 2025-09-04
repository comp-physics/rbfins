classdef TestCylinder < matlab.unittest.TestCase
    methods (Test)
        function smokeRun(testCase)
            % Set CI environment variable to ensure short run
            setenv('CI', 'true');
            
            try
                % Save current directory and add distmesh to path
                originalDir = pwd;
                addpath(fullfile(originalDir, 'distmesh'));
                
                % Test that the script runs without errors
                % We'll run it in eval to catch any syntax or runtime errors
                evalc('cylinder'); % evalc captures output and runs the script
                
                % If we get here, the script ran successfully
                fprintf('✅ Script executed successfully without errors!\n');
                testCase.verifyTrue(true, 'Cylinder script ran without errors');
                
            catch ME
                % If script fails, provide more informative error
                testCase.verifyTrue(false, sprintf('Cylinder script failed with error: %s', ME.message));
            end
        end
        
        function testConstants(testCase)
            % Test that all constants are properly defined in the script
            try
                % Read the script file and check for constant definitions
                scriptContent = fileread('cylinder.m');
                
                % Check that important constants are defined
                testCase.verifyTrue(contains(scriptContent, 'DOMAIN_X_MIN'), 'Domain X min constant should be defined');
                testCase.verifyTrue(contains(scriptContent, 'CYLINDER_RADIUS'), 'Cylinder radius constant should be defined');
                testCase.verifyTrue(contains(scriptContent, 'REYNOLDS_NUMBER'), 'Reynolds number constant should be defined');
                testCase.verifyTrue(contains(scriptContent, 'NUM_TIME_STEPS'), 'Number of time steps constant should be defined');
                testCase.verifyTrue(contains(scriptContent, 'STENCIL_SIZE_MAIN'), 'Main stencil size constant should be defined');
                
                fprintf('✅ All required constants are properly defined!\n');
                
            catch ME
                testCase.verifyTrue(false, sprintf('Constants test failed: %s', ME.message));
            end
        end
        
        function testNoMagicNumbers(testCase)
            % Test that magic numbers have been eliminated (spot check)
            try
                scriptContent = fileread('cylinder.m');
                
                % These magic numbers should no longer appear directly in the code
                % (we allow them in comments and constant definitions at the top)
                lines = strsplit(scriptContent, '\n');
                
                % Look for magic numbers in non-constant, non-comment lines
                magicNumberFound = false;
                problemLine = '';
                
                % Skip the first 100 lines (where constants are defined)
                for i = 101:length(lines)
                    line = lines{i};
                    % Skip comments and empty lines
                    trimmedLine = strtrim(line);
                    if isempty(trimmedLine) || startsWith(trimmedLine, '%')
                        continue;
                    end
                    
                    % Check for some specific magic numbers that should be replaced
                    if contains(line, 'k=35') || contains(line, 'k = 35')
                        magicNumberFound = true;
                        problemLine = line;
                        break;
                    end
                end
                
                testCase.verifyFalse(magicNumberFound, sprintf('Magic number found in line: %s', problemLine));
                
                fprintf('✅ Magic numbers have been properly replaced with constants!\n');
                
            catch ME
                testCase.verifyTrue(false, sprintf('Magic number test failed: %s', ME.message));
            end
        end
    end
end