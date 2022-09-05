close all
clear all
clc
% load("resultado_boderesistivo.mat")



sampling=0.1e-6;        % simulink periodo del paso de la simulacion en tiempo 
Duty=0.5;                    % simulink ciclo util del conversor
frec=logspace(1,4,20);         % simulink frecuencia de la señal pertubadora a la entrada
amp_bode=0.001;          % simulink amplitud de la señal perturbador a la entrada
sd_goal=0.01;               % desviación estandar maxima para decidir estado estable
db_goal=0.5;                  % amplitud de la perturbadora en la salida
time_vector=0:sampling:20e-3-sampling;
n=11;                          % numero de ciclos de la perturbadora a simular
%time_stop=((1/freq_bode)*n); % simulink tiempo de simulación de acuerdo a la señal perturbadora
time_start=0;
%max_data=time_stop/sampling% simulink maximo numero de datos a exportar a matlab
n_periodos_SS=10   % numero de periodos minimo para considerar estado estable;
margen=0.25;



freq_bode=1;
time_stop=round(((1/freq_bode)*n)-sampling,4); 
time_vector=0:sampling:time_stop;
amp_bode=0.002;
%%
out=sim('BodeConSimulink.slx','StartTime',num2str(time_start),'StopTime',num2str(time_stop)) ;

%%


%     for j=1:1:length(out.Signal_D.signals.values(1,1,:))
%         V_a2(j)=out.Signal_D.signals.values(1,1,j);
%     end
    %%
    hold on
    plot(V_a)
    plot(V_a2,"--")
%%

for k=1:1:length(frec)
    cont=1; 
    listo=0;
    time_stop=round(((1/frec(k))*n)-sampling,4); 
    ventana=1/(sampling*frec(k)); 
   %amp_bode=0.002;%(10.^((20*log10(0.014*(frec(k))))./20)+0.05)/1000;
    %n=resultado(k,5); %
   acumular=round((10.^((20*log10(100*(frec(k))))./20)+0.05)/1000);
  
    while (listo==0)
        time_vector=0:sampling:time_stop;
        time_stop=((1/frec(k))*n); 
        max_data=round(time_stop/sampling);
        freq_bode=frec(k);
        clear out
        out=sim('BodeConSimulink_perdidas.slx','StartTime',num2str(time_start),'StopTime',num2str(time_stop)) ;
        tam=length(out.Vcp_converter.signals.values);
        for i=1:1:n_periodos_SS
            promedios(i)=mean(out.Vcp_converter.signals.values((tam-round(ventana*i)):(tam-round(ventana*(i-1)))));
        end
        %plot(promedios)
        std_actual=std( promedios(end-n_periodos_SS+1:end) )
        
        % FFT
        V_a=out.Vcp_converter.signals.values(((tam-round(ventana*(n_periodos_SS-2))):(end)));   
        time_analisis=out.Vcp_converter.time(((tam-round(ventana*(n_periodos_SS-2))):(end)));
        T = sampling;             % Sampling period       
        L = length(time_analisis);             % Length of signal
        t = time_analisis;        % Time vector
        
        Y = fft(V_a);
        
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = (1/sampling)*(0:(L/2))/L;
        for j=1:1:length(f)
            if f(j)>=freq_bode
                pos=j;
                break;
            end
        end
        %Amp_funda=P1(pos);
        Amp_funda=max(V_a-mean(V_a))-min(V_a-mean(V_a));
        % plot(f,P1) 
        % title('Single-Sided Amplitude Spectrum of X(t)')
        % xlabel('f (Hz)')
        % ylabel('|P1(f)|')
        % xlim([0,10e3])
        if (std_actual<sd_goal)&&  (margen<Amp_funda) && (db_goal+margen>Amp_funda)
            listo=1;
        else 
            if std_actual>sd_goal
                n=n+acumular;
                %time_start=time_stop-(n_periodos_SS/frec(k))*5;
%                 if time_start<0
%                     time_start=0;
%                 end
            end
              if db_goal>Amp_funda
                amp_bode=amp_bode*3;
              end
               if db_goal+margen<Amp_funda
                amp_bode=amp_bode*0.5;
              end
    
        end
        if amp_bode<0.0001
            listo=1;
        end
        cont=cont+1;
        if cont==15
             listo=1;
        end
        amplitud_resultado(k)=amp_bode;
        disp(["Frecuencia: ",frec(k)]);
        disp(["Componente en frecuencia de la perturbadora: ",Amp_funda]);
        disp(["Desviación estandar del estado estable: ",std_actual]);
        disp(["Amplitud de la perturbadora: ",amp_bode]);
        %listo=1;
    end
   
    V_b=0.5+amp_bode*sin(2*pi*freq_bode*time_analisis);
    Amp_save(k,2)=Amp_funda;
    Amp_save(k,1)=amp_bode;
    Aavg=mean(V_a);
    Bavg=mean(V_b);
    Are=mean((V_a-Aavg).*cos(2*pi.*time_analisis.*freq_bode));
    Aim =mean((V_a-Aavg).*-sin(2*pi.*time_analisis.*freq_bode));
    Bre =mean((V_b-Bavg).*cos(2*pi.*time_analisis.*freq_bode));
    Bim =mean((V_b-Bavg).*-sin(2*pi.*time_analisis.*freq_bode));
    % corregir rsdianes y grados
    bode_result(k)=20*log10(hypot(Are,Aim) / hypot(Bre,Bim));
    phase_result(k)=mod(atan2(Aim, Are) - atan2(Bim, Bre)+pi,2*pi)-pi;
    Signal_In_A(k)=amp_bode;
    Signal_cont(k)=cont;
    signal_n(k)=n;
    plot(out.Vcp_converter.time,out.Vcp_converter.signals.values)

   % V_aTotal(:,k)=V_a;
   % V_bTotal(:,k)=V_b;
   if mod(i,5) ==0
    close all
    hold on 
    subplot(2,1,1)
    semilogx(frec(1:k),bode_result(1:k),"-o")
    subplot(2,1,2)
   semilogx(frec(1:k),phase_result(1:k),"-o")
   end
  if bode_result<0
      break;
  end
end
%%

%%
hold on
plot(V_b)
plot(V_a)
%%
hold on
plot((V_a-Aavg))
plot(cos(2*pi.*time_analisis.*freq_bode))

%%
close all
semilogx(frec,bode_result,"-o")
%%

close all
x=10:1:10e3;
y=20*log10(0.014*(x));
z=(10.^((y)./20)+0.05)/1000;
semilogx(x,z)
hold on
semilogx(frec(1:k-1),Signal_In_A(1:k-1),"-o")
legend("1","2")

%%

resultados=[frec;bode_result;phase_result]';
