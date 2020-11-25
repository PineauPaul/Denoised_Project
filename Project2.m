clear();

filename = "DDN1.wav";

%Recuperation of noiseless signal x
[x,f] = audioread("noiseless_samples/" + filename);
t_x = length(x)/f;

%Calculating number and length of signal fragments
t_x = length(x)/f;
nb_sig = floor(t_x/0.005); %4000
N = floor(length(x)/nb_sig); %330
divised = zeros(N,nb_sig);
x = x(1:N*nb_sig);


%Generation of noisy signal y
nb_noise = 100; %Number of fragments with noise only
x(1:nb_noise*N) = 0;
y = x;
y = awgn(y,10,'measured');
y_fft = fft(y);

%Division of y
for k = 1:nb_sig
    divised(:,k)=y(1+N*(k-1):N*k);
end 



%Step 1: estimate noise power spectral density
Sz = zeros(1,floor(N/2));
for k = 1:nb_noise
    z = divised(:,k);
    z_fft = fft(z);
    for i = 1:floor(N/2)
        Sz(i) = Sz(i) + (abs(z_fft(i))^2)/N;
    end   
end
Sz = Sz/nb_noise;



%Step 2: estimate noisy signal power spectral density
Sy = zeros(nb_sig,floor(N/2));
divised_fft = zeros(N,nb_sig);
for k = 1:nb_sig
    divised_fft(:,k) = fft(divised(:,k));
    for i = 1:floor(N/2)
        Sy(k,i) = (abs(divised_fft(i,k))^2)/N;
    end   
end


%Step 3: estimate noiseless signal power spectral density
Sx = zeros(nb_sig,floor(N/2));
for k = 1:nb_sig
   for i = 1:floor(N/2)
      if Sy(k,i)-Sz(i) > 0
          Sx(k,i) = Sy(k,i)-Sz(i);
      else
          Sx(k,i) = 0;
      end
   end
end


%Step 4: Design denoising filter
A = zeros(nb_sig,N);
for k = 1:nb_sig
   for i = 1:floor(N/2)
      A(k,i) = sqrt(Sx(k,i)/Sy(k,i));
   end
   for i = floor(N/2)+1:N
      A(k,i) = A(k,N+1-i);
   end
end



%Step 5: Design denoising filter
for k = 1:nb_sig
   denoised_frag_fft = zeros(N,1);
   for i = 1:N
      denoised_frag_fft(i) = A(k,i)*divised_fft(i,k);
   end
   denoised(1+N*(k-1):N*k) = real(ifft(denoised_frag_fft));
end


subplot(2, 3, 1);
plot(x);
title('x');
subplot(2, 3, 4);
plot(fft(x));
title('x fft');
subplot(2, 3, 2);
plot(y);
title('y');
subplot(2, 3, 5);
plot(fft(y));
title('y fft');
subplot(2, 3, 3);
plot(denoised);
title('denoised');
subplot(2, 3, 6);
plot(fft(denoised));
title('denoised fft');

audiowrite("noisy_signals\" + filename, y, f);
audiowrite("denoised_signals\" + filename, denoised, f);