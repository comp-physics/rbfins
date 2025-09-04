classdef TestCylinder < matlab.unittest.TestCase
    methods (Test)
        function smokeRun(testCase)
            % On GitHub Actions, CI=true is already set.
            % Run the main script; variables land in this workspace.
            cylinder

            % Basic existence checks
            testCase.verifyTrue(exist('W','var')==1);
            testCase.verifyTrue(exist('xy1','var')==1);

            % Dimensions consistent
            testCase.verifyEqual(size(W,1), 2*length(xy1));

            % No NaNs or Infs
            testCase.verifyFalse(any(isnan(W(:))));
            testCase.verifyFalse(any(isinf(W(:))));

            % Enforce no-slip on cylinder at final step (U,V zero on cylinder nodes)
            U = W(1:length(xy1), end);
            V = W(length(xy1)+1:end, end);
            testCase.verifyLessThan(max(abs(U(end-length(boundary_c)+1:end))), 1e-8);
            testCase.verifyLessThan(max(abs(V(end-length(boundary_c)+1:end))), 1e-8);
        end
    end
end
