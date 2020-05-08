%spatial step

Nr = 200; dr = 1/Nr; r = linspace(0,1,Nr+1).';

%time step 

dt = 1E-5; Nt = 100; %number of timesteps saved 

final_tau = 175; inner_time_steps = ceil(final_tau / (Nt*dt)); 

%create vectors for each variable

eta1 = ones(Nr+1, 1); eta2 = zeros(Nr+1,1); lambda = zeros(Nr+1, 1);

%dimensionless parameters via literature

delta1 = 12.5; rho2 = 1; delta2 = 4E-5; delta3 = 70;

%set initial conditions

eta1(end) = 0.01; eta2(end) = 1; lambda(end) = 1; %center of tumor

eta1(1) = 1; eta2(1) = 0; lambda(1) = 0;          %furthest from tumor

%create vectors for each time derivative 

eta1_dt = zeros(Nr+1,1); eta2_dt = zeros(Nr+1,1); lambda_dt = zeros(Nr+1,1);

%create vectors for each gradient and laplacian

eta1_dr = zeros(Nr-1,1); eta2_dr = zeros(Nr-1,1); lambda_dr = zeros(Nr-1,1);

eta2_dr2 = zeros(Nr-1,1); lambda_dr2 = zeros(Nr-1,1);

%tau matrix

tau = zeros(Nr+1,1);

%create matrix to save each value

eta1_matrix = zeros(Nt+1,Nr+1); eta2_matrix = zeros(Nt+1,Nr+1); 

lambda_matrix = zeros(Nt+1,Nr+1);

%documented time loop

for j = 1:Nt+1

%inner time loop

    for k = 1:inner_time_steps
    
        %compute gradients for each variable

        eta1_dr = compute_gradient(eta1,dr);

        eta2_dr = compute_gradient(eta2,dr);

        lambda_dr = compute_gradient(lambda,dr);

        %compute laplacian for eta2 and lambda

        eta2_dr2 = compute_laplacian(eta2,dr);

        lambda_dr2 = compute_laplacian(lambda,dr);

        % calculate time derivative for regular points

        eta1_dt(2:end-1) = eta1(2:end-1).*(1 - eta1(2:end-1)) -
            delta1*lambda(2:end-1).*eta1(2:end-1);

        eta2_dt(2:end-1) = rho2 * eta2(2:end-1).*(1 - eta2(2:end-1)) +
            delta2.*((1-eta1(2:end-1)).*((r(2:end-1).^4).*eta2_dr2(1:end) -eta2_dr(1:end).*
            (1-r(2:end-1)))-eta1_dr(1:end).*eta2_dr(1:end).*(r(2:end-1).^4));

        lambda_dt(2:end-1) = delta3*(eta2(2:end-1) - lambda(2:end-1)) +
            r(2:end-1).^3 .* (r(2:end-1).*lambda_dr2(1:end) +
            (2 - 1./(1-r(2:end-1))).* lambda_dr(1:end));

        %calculate time derivative for first point (r = 0)

        eta1_dt(end) = eta1(end) * (1 - eta1(end)) - delta1 * lambda(end) * eta1(end);
        
        eta2_dt(end) = rho2 * eta2(end) * (1 - eta2(end)) + 
            delta2 * (1 - eta1(end)) * r(end)^4 * (2 * eta2(end-1) - 2 * eta2(end)) / (dr^2);

        lambda_dt(end) = delta3 * (eta2(end) - lambda(end)) +
            r(end)^4 * (2 * lambda(end-1) - 2 * lambda(end)) / (dr^2);

        %calculate time derivative for last point (r=1)

        eta1_dt(1) = 0; eta2\_dt(1) = 0; lambda\_dt(1) = 0; 

        %increment time

        eta1 = eta1 + eta1_dt*dt; eta2 = eta2 + eta2_dt*dt;

        lambda = lambda + lambda_dt*dt;

    end

    %dose at tau = 25

     if (j == 25)

         D = 4; alpha = .8;

         eta1 = (eta1 + eta1_dt*dt).*exp(-alpha*D); 

         eta2 = (eta2 + eta2_dt*dt).* exp(-alpha*D);

     end
         
    %update matrices after inner_time_steps

    eta1_matrix(j,:) = eta1; eta2_matrix(j,:) = eta2;

    lambda_matrix(j,:) = lambda; tau(j) = dt*j*inner_time_steps;

end

save('output_chemo', 'eta1_matrix', 'eta2_matrix', 'lambda_matrix', 'tau','Nt', 'Nr', 'r'); \\