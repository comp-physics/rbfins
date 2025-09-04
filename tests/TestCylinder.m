classdef TestCylinder < matlab.unittest.TestCase
    methods (Test)
        function smokeRun(testCase)
            % This test simply verifies that the cylinder.m script can run
            % without errors when the CI environment variable is set
            
            % Set CI environment variable to ensure short run
            setenv('CI', 'true');
            
            % Set up paths
            setup_paths();
            
            % We'll use evalin('base') to run the script in the base workspace
            % This avoids the issue with the script clearing variables
            try
                evalin('base', 'success = false; errorMsg = ""; try, run(''cylinder.m''); success = true; catch ME, errorMsg = ME.message; end');
                
                % Check if the script ran successfully
                success = evalin('base', 'success');
                if success
                    fprintf('✅ Script executed successfully without errors!\n');
                    testCase.verifyTrue(true, 'Cylinder script ran without errors');
                else
                    errorMsg = evalin('base', 'errorMsg');
                    fprintf('❌ Script failed with error: %s\n', errorMsg);
                    testCase.verifyTrue(false, sprintf('Cylinder script failed with error: %s', errorMsg));
                end
            catch ME
                fprintf('❌ Test framework failed with error: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Test framework failed with error: %s', ME.message));
            end
        end
        
        function testConfigFile(testCase)
            % Test that the configuration file exists and contains all required parameters
            try
                % Set up paths
                setup_paths();
                
                % Check that config.m exists
                testCase.verifyTrue(exist('config.m', 'file') > 0, 'Configuration file config.m should exist');
                
                % Load the config file and check for required parameters
                cfg = config();
                
                % Check domain parameters
                testCase.verifyTrue(isfield(cfg, 'domain'), 'Config should have domain field');
                testCase.verifyTrue(isfield(cfg.domain, 'x_min'), 'Config should have domain.x_min');
                testCase.verifyTrue(isfield(cfg.domain, 'x_max'), 'Config should have domain.x_max');
                
                % Check mesh parameters
                testCase.verifyTrue(isfield(cfg, 'mesh'), 'Config should have mesh field');
                testCase.verifyTrue(isfield(cfg.mesh, 'cylinder_radius'), 'Config should have mesh.cylinder_radius');
                
                % Check simulation parameters
                testCase.verifyTrue(isfield(cfg, 'simulation'), 'Config should have simulation field');
                testCase.verifyTrue(isfield(cfg.simulation, 'reynolds_number'), 'Config should have simulation.reynolds_number');
                testCase.verifyTrue(isfield(cfg.simulation, 'num_time_steps'), 'Config should have simulation.num_time_steps');
                
                % Check RBF parameters
                testCase.verifyTrue(isfield(cfg, 'rbf'), 'Config should have rbf field');
                testCase.verifyTrue(isfield(cfg.rbf, 'stencil_size_main'), 'Config should have rbf.stencil_size_main');
                
                fprintf('✅ Configuration file contains all required parameters!\n');
                
            catch ME
                fprintf('❌ Config file test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Configuration file test failed: %s', ME.message));
            end
        end
        
        function testNoMagicNumbers(testCase)
            % Test that magic numbers have been eliminated (spot check)
            try
                % Set up paths
                setup_paths();
                
                scriptContent = fileread('cylinder.m');
                
                % These magic numbers should no longer appear directly in the code
                % Check for direct use of magic numbers that should be in config
                magicNumberPatterns = {
                    'k\s*=\s*35', 
                    'nu\s*=\s*1/100',
                    'dt\s*=\s*1e-2',
                    'Nt\s*=\s*5000',
                    'radius\s*=\s*0\.5'
                };
                
                % Look for magic numbers in non-constant, non-comment lines
                magicNumberFound = false;
                problemLine = '';
                
                lines = strsplit(scriptContent, '\n');
                for i = 1:length(lines)
                    line = lines{i};
                    % Skip comments and empty lines
                    trimmedLine = strtrim(line);
                    if isempty(trimmedLine) || startsWith(trimmedLine, '%')
                        continue;
                    end
                    
                    % Skip the constants section at the top
                    if contains(line, 'Configuration Constants') || i < 50
                        continue;
                    end
                    
                    % Check for magic number patterns
                    for p = 1:length(magicNumberPatterns)
                        if ~isempty(regexp(line, magicNumberPatterns{p}, 'once'))
                            magicNumberFound = true;
                            problemLine = line;
                            break;
                        end
                    end
                    
                    if magicNumberFound
                        break;
                    end
                end
                
                testCase.verifyFalse(magicNumberFound, sprintf('Magic number found in line: %s', problemLine));
                
                fprintf('✅ Magic numbers have been properly replaced with config parameters!\n');
                
            catch ME
                fprintf('❌ Magic number test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Magic number test failed: %s', ME.message));
            end
        end
        
        function testConfigUsage(testCase)
            % Test that cylinder.m uses the config file
            try
                % Set up paths
                setup_paths();
                
                scriptContent = fileread('cylinder.m');
                
                % Check that cylinder.m loads the config
                testCase.verifyTrue(contains(scriptContent, 'cfg = config()'), 'cylinder.m should load the config file');
                
                % Check that cylinder.m uses config parameters
                testCase.verifyTrue(contains(scriptContent, 'cfg.domain'), 'cylinder.m should use domain parameters from config');
                testCase.verifyTrue(contains(scriptContent, 'cfg.mesh'), 'cylinder.m should use mesh parameters from config');
                testCase.verifyTrue(contains(scriptContent, 'cfg.simulation'), 'cylinder.m should use simulation parameters from config');
                testCase.verifyTrue(contains(scriptContent, 'cfg.rbf'), 'cylinder.m should use RBF parameters from config');
                
                fprintf('✅ cylinder.m correctly uses the configuration file!\n');
                
            catch ME
                fprintf('❌ Config usage test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Config usage test failed: %s', ME.message));
            end
        end
    end
end