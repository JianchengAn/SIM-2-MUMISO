clc;
clearvars;
close all;
sigma2 = 10^(-104/10);
c = 3*10^8;
f0 = 28*10^9;
lambda = c/f0;
d_element_SIM = lambda/2;
MonteCarlo = 10;
n_max = 10;
N = 100;
K = 4;
Layer_Max = 10;
Pt = 10^(10/10)*10^(5/10);
R_MonteCarlo = zeros(MonteCarlo,1);
W_T = zeros(N, N);
Corr_T = zeros(N, N);
gama = zeros(K, 1);
R = zeros(K, 1);
location = [10*(1:K).' 10*(1:K).'];
W_T_1 = zeros(N, K);
R_layer = zeros(Layer_Max, 1);
gradient_temp = zeros(K, K);
del = zeros(K, 1);
for cc = 1 : Layer_Max
    L = cc;
    phase_transmit=zeros(N,L);
    d_layer = 5*lambda/L;
    for mm1 = 1 : N
        m_z = ceil(mm1/n_max);
        m_x = mod(mm1-1,n_max)+1;
        for mm2 = 1 : N
            n_z = ceil(mm2/n_max);
            n_x = mod(mm2-1,n_max)+1;
            d_temp  = sqrt( (m_x-n_x)^2 +  (m_z-n_z) ^2 )*d_element_SIM;
            d_temp2 = sqrt(d_layer^2 + d_temp^2);
            W_T(mm2,mm1) = lambda^2/4*(d_layer/d_temp2/d_temp2*(1/2/pi/d_temp2-1i/lambda))*exp(1i*2*pi*d_temp2/lambda);
            Corr_T(mm2,mm1) = sinc(2*d_temp/lambda);
        end
    end
    for mm = 1:N
        m_y = ceil(mm/n_max);
        m_x = mod(mm-1,n_max)+1;
        for nn = 1:K
            d_transmit = sqrt(d_layer^2 + ...
                ( (m_x-(1+n_max)/2)*d_element_SIM - (nn-(1+K)/2)*lambda/2)^2 + ...
                ( (m_y-(1+n_max)/2)*d_element_SIM )^2 );
            W_T_1(mm,nn) = lambda^2/4*(d_layer/d_transmit/d_transmit*(1/2/pi/d_transmit-1i/lambda))*exp(1i*2*pi*d_transmit/lambda);
        end
    end
    tic
    PD_transmit_phase = zeros(N,L);
    d = sqrt((10-5*lambda)^2 + location(:,1).^2+location(:,2).^2);
    pathloss = (lambda/(4*pi))^2./d.^(3.5);
    rng(1)
    for jj = 1:MonteCarlo
        G_in = sqrt(1/2)*(randn(K,N)+1i*randn(K,N));
        G = diag(sqrt(pathloss))*G_in*(Corr_T)^(1/2);
        phase_transmit = randn(N,L) + 1i*randn(N,L);
        phase_transmit = phase_transmit./abs(phase_transmit);
        W_SIM = diag(phase_transmit(:,1))*W_T_1;
        for l=1:L-1
            W_SIM = diag(phase_transmit(:,l+1))*W_T*W_SIM;
        end
        H_fit = G*W_SIM;
        p = IWF( Pt, sigma2, H_fit, K );
        for ii = 1:K
            gama(ii) = (abs(H_fit(ii, ii))^2*p(ii))/(abs(H_fit(ii, :)).^2*p - abs(H_fit(ii, ii))^2*p(ii) + sigma2);
            R(ii) = log2(1+gama(ii));
        end
        C_old = sum(R);
        C_new = C_old*2;
        phase_phase_transmit = angle(phase_transmit);
        count = 1;
        while abs(C_new-C_old) >= C_old * 0.000001 && count <= 50
            C_old = C_new;
            for ll = 1:L
                for mm = 1:N
                    X_left = W_T_1;
                    for ll_left = 1:ll-1
                        X_left = W_T*diag(phase_transmit(:,ll_left))*X_left;
                    end
                    X_right = G;
                    for ll_right = 1:(L-ll)
                        X_right = X_right*diag(phase_transmit(:,L+1-ll_right))*W_T;
                    end
                    for ss1 = 1:K
                        del(ss1) = 1/(abs(H_fit(ss1,:)).^2*p+sigma2);
                        for ss2 = 1:K
                            temp1 = X_right(ss1,mm)*X_left(mm,ss2);
                            if ss2 == ss1
                                gradient_temp(ss1,ss2) = del(ss1)*p(ss1)*imag((phase_transmit(mm,ll)*temp1)'*(H_fit(ss1,ss2)));
                            else
                                gradient_temp(ss1,ss2) = -del(ss1)*gama(ss1)*p(ss2)*imag((phase_transmit(mm,ll)*temp1)'*(H_fit(ss1,ss2)));
                            end
                        end
                    end
                    PD_transmit_phase(mm,ll) = 2/log(2)*sum(sum(gradient_temp));
                end
            end
            yy = pi/max(max(PD_transmit_phase));
            C_old_1 = C_old;
            C_new = 0;
            count_2 = 1;
            phase_transmit_temp = phase_transmit;
            while C_new < C_old_1 && count_2 <= 20
                phase_phase_transmit_temp = phase_phase_transmit+yy*PD_transmit_phase;
                phase_transmit_temp = exp(1i*phase_phase_transmit_temp);
                W_SIM = diag(phase_transmit_temp(:,1))*W_T_1;
                for l=1:L-1
                    W_SIM = diag(phase_transmit_temp(:,l+1))*W_T*W_SIM;
                end
                H_fit = G*W_SIM;
                p= IWF( Pt, sigma2, H_fit, K );
                for ii = 1:K
                    gama(ii) = (abs(H_fit(ii,ii))^2*p(ii))/(abs(H_fit(ii,:)).^2*p - abs(H_fit(ii,ii))^2*p(ii) + sigma2);
                    R(ii) = log2(1+gama(ii));
                end
                C_new = sum(R);
                yy = yy*0.5;
                count_2 = count_2+1;
            end
            phase_transmit = phase_transmit_temp;
            phase_phase_transmit = angle(phase_transmit);
            count = count+ 1;
        end
        R_MonteCarlo(jj) = C_new;
    end
    R_layer(cc) = mean(R_MonteCarlo);
    toc
end
figure
plot(R_layer)