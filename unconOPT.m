classdef unconOPT < handle

    properties
        %convergence parameter
        eps = 1e-6;
        x_opt
        f_opt
        i
    end

    methods

        function self = unconOPT(x,f,method)
            
            % reset iteration num
            self.i = 0;

            switch method
                case 'steepest-descend'
                    self.SteepestDescend(x,f);

                case 'classical-conjugate'
                    self.ClassicalConjugate(x,f);

                case 'hestenes-stiefel-conjugate'
                    self.HestenesStiefelConjugate(x,f);

                case 'fletcher-reeves-conjugate'
                    self.FletcherReevesConjugate(x,f);

                case 'polak-ribiere-conjugate'
                    self.PolakRibiereConjugate(x,f);

                case 'modified-newton'
                    self.ModifiedNewton(x,f);

                case 'davidon-fletcher-powell'
                    self.DFP(x,f);

                case 'broyden-fletcher-goldfarb-shanno'
                    self.BFGS(x,f);

            end
        end

        function SteepestDescend(self,x,f)

            % first gradient
            gradient = self.f_gradient(x,f);

            while norm(gradient) > self.eps
                d = -gradient;
                optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                if isnan(optStepSize)
                    return
                end
                x_new = x + optStepSize.*d;

                % calculating new gradient
                x = x_new;
                gradient = self.f_gradient(x,f);
                
                % store iteration solutions
                self.i = self.i + 1;
                self.x_opt(self.i,:) = x';
                self.f_opt(self.i) = f(x);
            end
        end


        function ClassicalConjugate(self,x,f)

            % first gradient
            gradient = self.f_gradient(x,f);
            first_iteration = true;

            while norm(gradient) > self.eps

                if first_iteration
                    d = -gradient;
                    optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                    if isnan(optStepSize)
                        return
                    end
                    x_new = x + optStepSize.*d;

                    % calculating new gradient
                    x = x_new;
                    new_gradient = self.f_gradient(x,f);
                    old_gradient = gradient;
                    gradient = new_gradient;

                    first_iteration = false;

                else
                    beta = norm(new_gradient)^2/norm(old_gradient)^2;
                    d = -new_gradient + beta*d;
                    optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                    if isnan(optStepSize)
                        return
                    end
                    x_new = x + optStepSize.*d;


                    x = x_new;
                    old_gradient = new_gradient;
                    new_gradient = self.f_gradient(x,f);
                    gradient = new_gradient;

                    % store iteration solutions
                    self.i = self.i + 1;
                    self.x_opt(self.i,:) = x';
                    self.f_opt(self.i) = f(x);
                end

            end
        end

        function HestenesStiefelConjugate(self,x,f)

            % first gradient
            gradient = self.f_gradient(x,f);
            first_iteration = true;

            while norm(gradient) > self.eps

                if first_iteration
                    d = -gradient;
                    optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                    if isnan(optStepSize)
                        return
                    end
                    x_new = x + optStepSize.*d;

                    % calculating new gradient
                    x = x_new;
                    new_gradient = self.f_gradient(x,f);
                    old_gradient = gradient;
                    gradient = new_gradient;

                    first_iteration = false;

                else
                    y = new_gradient - old_gradient;
                    beta = dot(new_gradient,y)/dot(d,y);
                    d = -new_gradient + beta*d;
                    optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                    if isnan(optStepSize)
                        return
                    end
                    x_new = x + optStepSize.*d;


                    x = x_new;
                    old_gradient = new_gradient;
                    new_gradient = self.f_gradient(x,f);
                    gradient = new_gradient;

                    % store iteration solutions
                    self.i = self.i + 1;
                    self.x_opt(self.i,:) = x';
                    self.f_opt(self.i) = f(x);

                end

            end

        end


        function FletcherReevesConjugate(self,x,f)
            self.ClassicalConjugate(x,f)
        end

        function PolakRibiereConjugate(self,x,f)

            % first gradient
            gradient = self.f_gradient(x,f);
            first_iteration = true;

            while norm(gradient) > self.eps

                if first_iteration
                    d = -gradient;
                    optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                    if isnan(optStepSize)
                        return
                    end
                    
                    x_new = x + optStepSize.*d;

                    % calculating new gradient
                    x = x_new;
                    new_gradient = self.f_gradient(x,f);
                    old_gradient = gradient;
                    gradient = new_gradient;

                    first_iteration = false;

                else
                    y = new_gradient - old_gradient;
                    beta = dot(new_gradient,y)/dot(old_gradient,old_gradient);
                    d = -new_gradient + beta*d;
                    optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                    if isnan(optStepSize)
                        return
                    end
                    x_new = x + optStepSize.*d;

                    x = x_new;
                    old_gradient = new_gradient;
                    new_gradient = self.f_gradient(x,f);
                    gradient = new_gradient;

                    % store iteration solutions
                    self.i = self.i + 1;
                    self.x_opt(self.i,:) = x';
                    self.f_opt(self.i) = f(x);
                end
            end

        end


        function ModifiedNewton(self,x,f)

            gradient = self.f_gradient(x,f);

            while norm(gradient) > self.eps

                hessian = self.f_hessian(x,f);
                eig_hessian = eig(hessian);

                % checking Hessian is semi definite
                % if ~all(eig_hessian > 0)
                %     disp('Search is terminated. Because Hessian is not semi-positive definite')
                %     self.x_opt = x;
                %     return
                % end

                d = linsolve(hessian,-gradient);
                optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                if isnan(optStepSize)
                    return
                end
                x_new = x + optStepSize.*d;

                % calculating new gradient
                x = x_new;
                gradient = self.f_gradient(x,f);


                % store iteration solutions
                self.i = self.i + 1;
                self.x_opt(self.i,:) = x';
                self.f_opt(self.i) = f(x);
            end
        end

        function DFP(self,x,f)

            A = eye(length(x));

            gradient = self.f_gradient(x,f);

            while norm(gradient) > self.eps

                d = -A*gradient;
                optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                if isnan(optStepSize)
                    return
                end
                x_new = x + optStepSize.*d;
                x = x_new;
                new_gradient = self.f_gradient(x,f);


                % update A
                s = optStepSize * d;
                y = new_gradient - gradient;
                z = A * y;
                B = (s*s') / (dot(s,y));
                C = (-z*z') / (dot(y,z));
                A = A + B + C;

                % update gradient
                gradient = new_gradient;

                % store iteration solutions
                self.i = self.i + 1;
                self.x_opt(self.i,:) = x';
                self.f_opt(self.i) = f(x);
            end
        end

        function BFGS(self,x,f)

            H = eye(length(x));

            gradient = self.f_gradient(x,f);

            while norm(gradient) > self.eps

                d = linsolve(H,-gradient);
                optStepSize = findOptStepSizeGoldenSearch(f, x, d);
                x_new = x + optStepSize.*d;
                x = x_new;
                new_gradient = self.f_gradient(x,f);


                % update H
                s = optStepSize * d;
                y = new_gradient - gradient;
                D = (y*y') / (dot(y,s));
                E = (gradient*gradient') / (dot(gradient,d));
                H = H + D + E;

                % update gradient
                gradient = new_gradient;

                % store iteration solutions
                self.i = self.i + 1;
                self.x_opt(self.i,:) = x';
                self.f_opt(self.i) = f(x);
            end
        end
    end

    methods(Static)

        function gradient = f_gradient(x,f)

            eps = 1e-12;
            x1 = x;
            x2 = x;

            gradient = zeros(length(x1),1);

            for i=1:length(x1)
                x2(i) = x(i) + eps;
                gradient(i) = (f(x2) - f(x1) ) / norm(x2 - x1);
                x2 = x;
            end
        end

        function hessian = f_hessian(x,f)

            n = length(x);
            hessian = zeros(n, n);
            eps = 1e-6; % small value for numerical differentiation

            % Compute the diagonal elements of the Hessian matrix
            % for i = 1:n
            %     dx = zeros(n, 1);
            %     dx(i) = eps;
            %     f1 = f(x + dx);
            %     f2 = f(x - dx);
            %     H(i, i) = (f1 - 2*f(x) + f2) / (eps^2);
            % end

            % Compute the off-diagonal elements of the Hessian matrix
            for i = 1:n
                for j = 1:n
                    dx1 = zeros(n, 1);
                    dx2 = zeros(n, 1);
                    dx1(i) = eps;
                    dx2(j) = eps;
                    f1 = f(x + dx1 + dx2);
                    f2 = f(x + dx1 - dx2);
                    f3 = f(x - dx1 + dx2);
                    f4 = f(x - dx1 - dx2);
                    hessian(i, j) = (f1 - f2 - f3 + f4) / (4 * eps^2);
                    hessian(j, i) = hessian(i, j);
                end
            end
        end

    end
end
