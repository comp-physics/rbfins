classdef TestCylinder < matlab.unittest.TestCase
    methods (Test)
        function smokeRun(testCase)
            % Test that the simulate.m script runs without errors in CI environment
            
            % Set CI environment variable to ensure short run
            setenv('CI', 'true');
            setenv('MATLAB_TEST', 'true');
            
            % Set up paths
            setup_paths();
            
            % Test that the script runs successfully
            try
                % Run the script and capture any errors
                evalin('base', 'clear; success = false; errorMsg = "";');
                evalin('base', 'try; simulate; success = true; catch ME; errorMsg = ME.message; end;');
                
                % Check if the script ran successfully
                success = evalin('base', 'success');
                
                if success
                    fprintf('✅ Cylinder script executed successfully without errors!\n');
                    testCase.verifyTrue(true, 'Cylinder script completed successfully');
                else
                    errorMsg = evalin('base', 'errorMsg');
                    fprintf('❌ Cylinder script failed with error: %s\n', errorMsg);
                    testCase.verifyTrue(false, sprintf('Cylinder script failed: %s', errorMsg));
                end
                
            catch ME
                fprintf('❌ Test execution failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Test execution failed: %s', ME.message));
            end
        end
        
        function testConfigFileExists(testCase)
            % Test that the configuration file exists and loads correctly
            try
                setup_paths();
                
                % Check that config.m exists
                testCase.verifyTrue(exist('config.m', 'file') > 0, 'config.m file should exist');
                
                % Test that config function can be called
                cfg = config();
                testCase.verifyTrue(isstruct(cfg), 'config() should return a structure');
                
                % Check for essential fields
                testCase.verifyTrue(isfield(cfg, 'domain'), 'Config should have domain field');
                testCase.verifyTrue(isfield(cfg, 'mesh'), 'Config should have mesh field');
                testCase.verifyTrue(isfield(cfg, 'simulation'), 'Config should have simulation field');
                
                fprintf('✅ Configuration file test passed!\n');
                
            catch ME
                fprintf('❌ Configuration file test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Configuration file test failed: %s', ME.message));
            end
        end
        
        function testMainScriptStructure(testCase)
            % Test that the main script has the expected structure
            try
                setup_paths();
                
                % Check that simulate.m exists
                testCase.verifyTrue(exist('simulate.m', 'file') > 0, 'simulate.m file should exist');
                
                % Read the script content
                content = fileread('simulate.m');
                
                % Check for key components
                testCase.verifyTrue(contains(content, 'cfg = config()'), 'Script should load configuration');
                testCase.verifyTrue(contains(content, 'isCI'), 'Script should check CI environment');
                testCase.verifyTrue(contains(content, 'addpath'), 'Script should add distmesh path');
                
                fprintf('✅ Main script structure test passed!\n');
                
            catch ME
                fprintf('❌ Main script structure test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Main script structure test failed: %s', ME.message));
            end
        end
        
        function testDependenciesExist(testCase)
            % Test that required dependencies exist
            try
                setup_paths();
                
                % Check for essential files
                testCase.verifyTrue(exist('NS_2d_cylinder_PHS.m', 'file') > 0, 'NS_2d_cylinder_PHS.m should exist');
                testCase.verifyTrue(exist('RBF_PHS_FD_all.m', 'file') > 0, 'RBF_PHS_FD_all.m should exist');
                testCase.verifyTrue(exist('nearest_interp.m', 'file') > 0, 'nearest_interp.m should exist');
                testCase.verifyTrue(exist('knnsearch.m', 'file') > 0, 'knnsearch.m should exist');
                
                % Check for distmesh directory
                testCase.verifyTrue(exist('distmesh', 'dir') > 0, 'distmesh directory should exist');
                
                fprintf('✅ Dependencies test passed!\n');
                
            catch ME
                fprintf('❌ Dependencies test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Dependencies test failed: %s', ME.message));
            end
        end
    end
end