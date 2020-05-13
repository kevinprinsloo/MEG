
cd /analyse/Project0129

%% localise evoked component in all 8 conditions to identify ROIs
allS=zeros(14,8,5798,408);
for is=1:14,
    isub=goodsub(is)
        for icond=1:8,
            cd ../cleanDataICA
                loadname=['datacl2_' num2str(isub) '_' num2str(icond) '.mat'];
                load(loadname)
               load(['../GradInfo/data_' num2str(isub) '_' num2str(icond) '.mat'])
                 data3.grad=grad;
                cfg=[];
                if icond<5,
                cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
                else
                   cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan13{:}'])};
                end  
                data3=ft_preprocessing(cfg,data3);
                NT=length(data3.trial);
                for k=1:NT,
                    data3.trial{k}=data3.trial{k}*1e13;
                end
           
            
            cfg=[];
            cfg.lpfilter='yes';
            cfg.lpfreq=40;
            %cfg.hilbert='complex';
            tmp=ft_preprocessing(cfg,data3);
            
              cfg=[];
            cfg.toilim=[-0.1 0.3];
            databp=ft_redefinetrial(cfg,tmp);
            
            cfg=[];
            cfg.covariance='yes';
            cfg.covariancewindow=[-0.1 0.3];
            avg=ft_timelockanalysis(cfg,databp);

          

            cd ../Source
            %outname=['MNE_filter' num2str(is) '-cond' num2str(3)];
            %load(outname)
            
            
            load(['Vol_' num2str(isub)])
            load(['Grid_' num2str(isub)])
             
            cfg=[];
            cfg.method= 'eloreta';
            cfg.eloreta.keepfilter='yes';
            cfg.eloreta.reducerank='yes';
            cfg.eloreta.normalize='yes';
            cfg.grid= indgrid;
            cfg.vol= vol;
            if icond<5,
                cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
                else
                   cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan13{:}'])};
                end  
            cfg.grad=avg.grad;
            cfg.reducerank=2;
            cfg.normalize='yes';
            grid=ft_prepare_leadfield(cfg,avg);
            cfg.grid=grid;
            source = ft_sourceanalysis(cfg, avg);
            cd ../Kevin
            filt=cat(1,source.avg.filter{source.inside});
            ori=cat(1,source.avg.ori{source.inside});
            %po=cat(1,source.avg.mom{source.inside});
            %allS(is,icond,:,:)=po(1:2:end,:).^2+po(2:2:end,:).^2;
            outname=['EL_filter' num2str(is) '-cond' num2str(icond)];
            save(outname,'filt','ori')
        end
end

allS2=repmat(mean(allS(:,:,:,1:40),4),[1 1 1 408]);  %use mean values over time

tmp=squeeze(allS-allS2);
me=squeeze(mean(tmp,1));
st=squeeze(std(tmp,[],1))/sqrt(nsubj);
tst=me./st;clear tmp me st allS2

tmp=squeeze(allS(:,2,:,:)-allS(:,4,:,:));
me=squeeze(mean(tmp,1));
st=squeeze(std(tmp,[],1))/sqrt(nsubj);
tst=me./st;clear tmp me st allS2


%% localise phase resetting
cd /analyse/Project0129
subinfor

cd Source


allS=zeros(4,5,14,5798,102);
ff=[1 3;3 7; 7 13;15 25; 30 45];
goodsub=[1:14];
for is=1:14,
    isub=goodsub(is)
    
    for icond=1:4,
        for f1=1:5,
            cd ../cleanDataICA
            if (f1==1),
                loadname=['datacl2_' num2str(isub) '_' num2str(icond) '.mat'];
                load(loadname)
           load(['../GradInfo/data_' num2str(isub) '_' num2str(icond) '.mat'])
                 data3.grad=grad;
                cfg=[];
                cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
                data3=ft_preprocessing(cfg,data3);
                NT=length(data3.trial);
                for k=1:NT,
                    data3.trial{k}=data3.trial{k}*1e13;
                end
            end
            
            cfg=[];
            cfg.bpfilter='yes';
            if (f1 ==1)
                cfg.bpfiltord=3;
            end
            cfg.bpfreq=[ff(f1,1) ff(f1,2)];
            %cfg.hilbert='complex';
            tmp=ft_preprocessing(cfg,data3);
            
            
            cfg=[];
            cfg.covariance='yes';
            cfg.covariancewindow=[-0.5 0.5];
            avg=ft_timelockanalysis(cfg,tmp);
            
            cfg=[];
            cfg.hilbert='complex';
            databp=ft_preprocessing(cfg,tmp);
            cfg=[];
            cfg.toilim=[-0.5 0.5];
            databp=ft_redefinetrial(cfg,databp);

            cd ../Source
            %outname=['MNE_filter' num2str(is) '-cond' num2str(3)];
            %load(outname)
            
            
            load(['Vol_' num2str(isub)])
            load(['Grid_' num2str(isub)])
             
            cfg=[];
            cfg.method= 'lcmv';
            cfg.lcmv.lambda='5%';
            cfg.lcmv.keepfilter='yes';
            cfg.lcmv.reducerank='yes';
            cfg.lcmv.normalize='yes';
            cfg.grid= indgrid;
            cfg.vol= vol;
            cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
            cfg.grad=avg.grad;
            cfg.reducerank=2;
            cfg.normalize='yes';
            if (f1==1),
                grid=ft_prepare_leadfield(cfg,avg);
            end;
            cfg.grid=grid;
            source = ft_sourceanalysis(cfg, avg);
            filt=cat(1,source.avg.filter{source.inside})*1e12;
            
            NT=length(databp.trial);
            %
            NS=length(databp.time{1});
            PH=zeros(3,NT,5798,102);
                         for k=1:NT,
                            for k2=1:3,
                                PH(k2,k,:,:)=angle(filt(k2:3:end,:)*databp.trial{k}(:,1:10:end));
                            end
                         end
          
            
            %for EL only two filters
            % for k=1:NT,
            %   t1=abs(filt(1:2:end,:)*databp.trial{k}(:,1:10:end)).^2+abs(filt(2:2:end,:)*databp.trial{k}(:,1:10:end)).^2;

                %t1=abs(filt(1:3:end,:)*databp.trial{k}).^2+abs(filt(2:3:end,:)*databp.trial{k}).^2+abs(filt(3:3:end,:)*databp.trial{k}).^2;
               % t1=t1*1e30;
             %   bb=repmat(mean(t1,2),[1 101]);
              %  t2=t2+(t1-bb)./bb;
                %allS(icond,f1,is,:,:)=squeeze(allS(icond,f1,is,:,:))+10*log10(t1./repmat(mean(t1,2),[1 102]));
             allS(icond,f1,is,:,:)=mean(abs(mean(exp(i*PH),2)),1);

        end
           
        end
end

%visual?
figure
plot(ti,squeeze(mean(allS(:,3,:,1548,:),3)));legend({'1' '2' '3' '4'})
figure
plot(ti,squeeze(mean(allS(:,3,:,1906,:),3)))
legend({'1' '2' '3' '4'})
%right auditory
figure
plot(ti,squeeze(mean(allS(:,3,:,2916,:),3)))
legend({'1' '2' '3' '4'})

allS2=repmat(mean(allS(:,:,:,:,1:40),5),[1 1 1 1 102]);  %use mean values over time

tmp=squeeze(allS-allS2);
me=squeeze(mean(tmp,3));
st=squeeze(std(tmp,[],3))/sqrt(nsubj);
tst=me./st;clear tmp me st allS2




goodsub=[1:14];nsubj=length(goodsub);
tmp=squeeze((allS(1,:,goodsub,:,:))-(allS(3,:,goodsub,:,:)));
me=squeeze(mean(tmp,2));
st=squeeze(std(tmp,[],2))/sqrt(nsubj);
tst=me./st;clear tmp me st

load /data1/synchro1/Software/fieldtrip-20131217/template/sourcemodel/standard_sourcemodel3d8mm.mat
mri=ft_read_mri('/data1/synchro1/Software/fieldtrip-20131217/template/anatomy/single_subj_T1_1mm.nii');

source2=source;
source2=rmfield(source2,'time');
source2.avg.pow=zeros(11000,1);
source2.pos=sourcemodel.pos;

for k=60:2:100,
    source2.avg.pow(source2.inside)=squeeze(tst(4,3,:,k));
 cfg=[];
 cfg.downsample=2;
 cfg.parameter={'avg.pow'};
 interp=ft_sourceinterpolate(cfg,source2,mri);
 close all
 cfg=[];
 cfg.funparameter='avg.pow';
  cfg.method='slice';
 cfg.interactive='yes';
 %cfg.funcolorlim=[-2 2];
 ft_sourceplot(cfg,interp);
 title(k)
 pause
 
end
 
figure
plot(ti,squeeze(mean(allS(:,3,:,1906,:),3)))
legend({'1' '2' '3' '4'})

%% stats on ITC with module diff_ITC?







%% look at TFR of phase and power of few voxels
ro=[5740 5051 6618 7094]; %aud Cort, first three right
ro=[5757 4758 4739 6224 5224 4743 8206 9708 8197 9715 6113 4571 5572 4551 5532 1168 1175];
rolab={'rAC', 'rAC', 'rAC', 'lAC' 'lAC' 'lAC', 'LM1' 'LMI' 'RM1' 'RM1', 'vis', 'vis2', 'vis3', 'vis4', 'vis5', 'cer1' 'cer2' 'cer3'};

%from stats
%Vis 6113 4571 5572 4551 5532
 nro=length(ro);
for k=1:nro,
    iro(k)=find(source.inside == ro(k));
end

%allS=zeros(4,5,14,nro,2,2,101);
goodsub=[1:14];
allPO=zeros(14,4,nro,50,178);
allPL=allPO;
for is=1:14,
    is
    isub=goodsub(is)
    
    for icond=1:4,
        cd ../cleanDataICA
        loadname=['datacl2_' num2str(isub) '_' num2str(icond) '.mat'];
        load(loadname)
        
        cfg=[];
        if icond<5,
                cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
                else
                   cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan13{:}'])};
                end  
        data3=ft_preprocessing(cfg,data3);
        
        cfg=[];
        cfg.toilim=[-1 1];
        data3=ft_redefinetrial(cfg,data3);
        
        cd ../Kevin
        outname=['EL_filter' num2str(is) '-cond' num2str(icond)];
        load(outname)
        ori2=reshape(ori,2,11596/2);
        
        NT=length(data3.trial);
        %
        NS=length(data3.time{1});
        
%         dat=zeros(nro,NT,2,2035);
%         for ir=1:nro,
%             for k=1:NT
%                 for k2=1:2,
%                     dat(ir,k,k2,:)=filt(iro(ir)*2-2+k2,:)*(data3.trial{k}*1e13);
%                 end
%             end
%         end
        
%use optimal ori
          dat=zeros(nro,NT,2035);
        for ir=1:nro,
            for k=1:NT
              tmp1=filt(iro(ir)*2-1,:)*(data3.trial{k}*1e13);
              tmp2=filt(iro(ir)*2,:)*(data3.trial{k}*1e13);
             
                dat(ir,k,:)=ori2(iro(ir)*2-1)*tmp1+ori2(iro(ir)*2)*tmp2;
               
            end
        end
        
        %do tfr
        for ir=1:nro,
            %keep single trial tfr
            tt=zeros(NT,50,178);
            for k=1:NT,
                %for oo=1:2,
                    %tmp=ori{ir}(:,oo)'*squeeze(dat(ir,k,:,:));
                    tmp=squeeze(dat(ir,k,:));
                    %[y,f,t] = spectrogram(tmp,256,250,[1:40],1017); % uses DFT
                    [y,f,t] = spectrogram(tmp,256,246,512,1017); % uses FFT
                    %cw=cwtft({tmp,1/1017},'scales',{0.0097, 0.12,50,'pow'});
                    tt(k,:,:)=y(1:50,:);
                    %tt(k,:,:)=cw.cfs(:,550:5:1435);
                %end
            end
            %compute power and plv
            t=squeeze(abs(tt).^2);
            %t2=squeeze(mean(((t-repmat(mean(t,3),[1 1 178]))./repmat(mean(t,3),[1 1 178])),1));
            t2=squeeze(mean(10*log10(t./repmat(mean(t,3),[1 1 178])),1));

            %t2=mean(t,1);
            allPO(is,icond,ir,:,:)=squeeze(t2);
            allPL(is,icond,ir,:,:)=squeeze(abs(mean(exp(i*angle(tt)),1)));
%             for k3=1:50,
%                 for k4=1:178,
%                     allPL(is,icond,ir,k3,k4)=circ_mean(squeeze(angle(tt(1,:,k3,k4))));
%                 end
%             end
        end
    end
end



%% significant phase difference aud cortex and vis cortex
allPO=zeros(14,4,nro,50,178);
allPL=allPO;
%combinations lf AC voxel and vis voxels to compute
m=1;
for k1=1:6,
    for k2=11:15,
cv(m,1)=k1;
cv(m,2)=k2;
    m=m+1;
    end;
end

for ic=24:30,
for is=1:14,
    
    isub=goodsub(is)
    
    for icond=4:4,
        cd ../cleanDataICA
        loadname=['datacl2_' num2str(isub) '_' num2str(icond) '.mat'];
        load(loadname)
        
        cfg=[];
        if icond<5,
                cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan' num2str(icond) '{:}'])};
                else
                   cfg.channel={'MEG',eval(['SUBJ{' num2str(isub) '}.badchan13{:}'])};
                end  
        data3=ft_preprocessing(cfg,data3);
        
        cfg=[];
        cfg.toilim=[-1 1];
        data3=ft_redefinetrial(cfg,data3);
        
        cd ../Kevin
        outname=['EL_filter' num2str(is) '-cond' num2str(icond)];
        load(outname)
        ori2=reshape(ori,2,11596/2);
        
        NT=length(data3.trial);
        %
        NS=length(data3.time{1});
        
%         dat=zeros(nro,NT,2,2035);
%         for ir=1:nro,
%             for k=1:NT
%                 for k2=1:2,
%                     dat(ir,k,k2,:)=filt(iro(ir)*2-2+k2,:)*(data3.trial{k}*1e13);
%                 end
%             end
%         end
        
%use optimal ori
          dat=zeros(nro,NT,2035);
        for ir=1:nro,
            for k=1:NT
              tmp1=filt(iro(ir)*2-1,:)*(data3.trial{k}*1e13);
              tmp2=filt(iro(ir)*2,:)*(data3.trial{k}*1e13);
             
                dat(ir,k,:)=ori2(iro(ir)*2-1)*tmp1+ori2(iro(ir)*2)*tmp2;
               
            end
        end
        
        %do tfr
            %keep single trial tfr
            tt=zeros(NT,50,178);
            tt2=tt;
            for k=1:NT,
                %for oo=1:2,
                    %tmp=ori{ir}(:,oo)'*squeeze(dat(ir,k,:,:));
                    tmp=squeeze(dat(cv(ic,1),k,:));
                    %[y,f,t] = spectrogram(tmp,256,250,[1:40],1017); % uses DFT
                    [y,f,t] = spectrogram(tmp,256,246,512,1017); % uses FFT
                    %cw=cwtft({tmp,1/1017},'scales',{0.0097, 0.12,50,'pow'});
                    tt(k,:,:)=y(1:50,:);
                    tmp=squeeze(dat(cv(ic,2),k,:));
                    %[y,f,t] = spectrogram(tmp,256,250,[1:40],1017); % uses DFT
                    [y,f,t] = spectrogram(tmp,256,246,512,1017); % uses FFT
                    %cw=cwtft({tmp,1/1017},'scales',{0.0097, 0.12,50,'pow'});
                    tt2(k,:,:)=y(1:50,:);
                    
                    %tt(k,:,:)=cw.cfs(:,550:5:1435);
                %end
            end
            allNT(is)=NT;
            if (is==1)
                alltt=angle(tt);
                alltt2=angle(tt2);
            else
                alltt=[alltt;angle(tt)];
                alltt2=[alltt2;angle(tt2)];
                
            end
            %compute power and plv
%             t=squeeze(abs(tt).^2);
%             %t2=squeeze(mean(((t-repmat(mean(t,3),[1 1 178]))./repmat(mean(t,3),[1 1 178])),1));
%             t2=squeeze(mean(10*log10(t./repmat(mean(t,3),[1 1 178])),1));
% 
%             %t2=mean(t,1);
%             allPO(is,icond,ir,:,:)=squeeze(t2);
%             allPL(is,icond,ir,:,:)=squeeze(abs(mean(exp(i*angle(tt)),1)));
%             for k3=1:50,
%                 for k4=1:178,
%                     allPL(is,icond,ir,k3,k4)=circ_mean(squeeze(angle(tt(1,:,k3,k4))));
%                 end
%             end
        end
    end
  

    %cmtest
%r=zeros(30,50,178);
pp=zeros(50,178);
for k1=1:20,
    for k2=80:160,
%         [r(1,k1,k2),r(2,k1,k2),r(3,k1,k2)]=circ_mean(squeeze(alltt(:,k1,k2))');
%         if (r(2,k1,k2) <0 || r(3,k1,k2)>0)
%         pp(k1,k2)=1;
%         end
r(ic,k1,k2)=circ_cmtest(squeeze(alltt(:,k1,k2))',squeeze(alltt2(:,k1,k2))');
    end
end
end

%saved as CM_test2  /analyse/Project0129/Kevin
H=fspecial('gaussian',[2 5],2);
for k=1:30,
tmp=imfilter(squeeze(r(k,14,:,:)),H,'replicate');

for k=[7 26],
    figure(1)
    imagesc(t-1,f(1:50),squeeze(mean(-log10(r(k,:,:)),1)));colorbar;
    title(k);
    figure(2)
    tmp=imfilter(squeeze(r(k,:,:)),H,'replicate');
    imagesc(t-1,f(1:50),squeeze(-log10(tmp)));colorbar;

    pause;
end;

  


%for each individual
r=zeros(14,50,178);
m=1;
for s=1:14,
    s
    ti=[m:m+allNT(s)];
    m=m+allNT(s);
for k1=1:20,
    for k2=80:178,
%         [r(1,k1,k2),r(2,k1,k2),r(3,k1,k2)]=circ_mean(squeeze(alltt(:,k1,k2))');
%         if (r(2,k1,k2) <0 || r(3,k1,k2)>0)
%         pp(k1,k2)=1;
%         end
r(s,k1,k2)=circ_cmtest(squeeze(alltt(ti,k1,k2))',squeeze(alltt2(ti,k1,k2))');
    end
end
end

%%
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
freq=1./(cw.scales*MorletFourierFactor);

nsubj=14;

allS2=repmat(mean(allPL(:,:,:,:,1:80),5),[1 1 1 1 178]);  %use mean values over time

tmp=squeeze(allPL-allS2);
me=squeeze(mean(tmp,1));
st=squeeze(std(tmp,[],1))/sqrt(nsubj);
tst=me./st;clear tmp me st allS2


%adjust for overall power difference between blocks
%allPO=allPO-repmat(mean(allPO(:,:,:,:,:,:),6),[1 1 1 1 1 178]);
tmp=squeeze(mean(allPO(goodsub,1,:,:,:)-allPO(goodsub,3,:,:,:),2));  %maybe log10
me=squeeze(mean(tmp,1));
st=squeeze(std(tmp,[],1))/sqrt(nsubj);
tst=me./st;clear tmp me st

tmp=squeeze(mean(allPL(:,1,:,:,:)-allPL(:,3,:,:,:),2));
me=squeeze(mean(tmp,1));
st=squeeze(std(tmp,[],1))/sqrt(nsubj);
tst2=me./st;clear tmp me st

t=t-1;

for k=1:14,
    pcolor(tim,f(1:50),squeeze(mean(allPO(k,1,3,:,:),1)));shading interp;colorbar
    title(k)
    pause
end

for k=11:15,
    pcolor(tim,f(1:50),squeeze(mean(tst(4,k,:,:),1)));shading interp;colorbar
    title(rolab{k})
    pause
end


pcolor(tim,f(1:50),squeeze(mean(tst(4,14,:,:),1)));shading interp;colorbar
tmp=imfilter(squeeze(tst(4,14,:,:)),H,'replicate');

%check t-test agains aud
for k=1:6,
tmp=squeeze(mean(allPL(:,4,14,:,:)-allPL(:,4,k,:,:),2));
me=squeeze(mean(tmp,1));
st=squeeze(std(tmp,[],1))/sqrt(nsubj);
tst2=me./st;clear tmp me st
pcolor(tim,f(1:50),squeeze((tst2(:,:))));shading interp;colorbar
pause;
end



for k=1:15,plot(tim,squeeze(mean(allPL(:,1:4,k,6,:),1)));title(k);pause;end;
%correlation
co=zeros(15,50,178);
co2=co;
for k1=1:15,
    k1
    for k2=1:50,
        for k3=1:178,
            co(k1,k2,k3)=corr(atten,squeeze(allPL(:,1,k1,k2,k3)-allPL(:,3,k1,k2,k3)));
            co2(k1,k2,k3)=corr(atten,squeeze(allPO(:,1,k1,k2,k3)-allPO(:,3,k1,k2,k3)));
            
        end
    end
end


%% TFR on ICA components
tmp=squeeze(comp.trial{1312+1}(1,:));
w1=60;ov=w1-3;fftsize=256;
%[y,f,t] = spectrogram(tmp,w1,ov,[1:40],250); % uses DFT
[y,f,t] = spectrogram(tmp,w1,ov,fftsize,250); % uses FFT

allPO=zeros(30,30,size(y,1),size(y,2));
allPOa=zeros(30,size(y,1));
allPL=allPO;
for ic=1:30, %1312
    ic
     tt=zeros(1315,size(y,1),size(y,2));
    for k=1:1315,
        
        %tmp=ori{ir}(:,oo)'*squeeze(dat(ir,k,:,:));
        tmp=squeeze(comp.trial{1312+k}(ic,:));
        %[y,f,t] = spectrogram(tmp,w1,ov,[1:40],250); % uses DFT
        [y,f,t] = spectrogram(tmp,w1,ov,fftsize,250); % uses FFT
        tt(k,:,:)=y;
    end
    for ic2=ic+1:30,
    %keep single trial tfr
    tt2=zeros(1315,size(y,1),size(y,2));
    for k=1:1315,
        
        %tmp=ori{ir}(:,oo)'*squeeze(dat(ir,k,:,:));
        tmp=squeeze(comp.trial{1312+k}(ic2,:));
        %[y,f,t] = spectrogram(tmp,w1,ov,[1:40],250); % uses DFT
        [y,f,t] = spectrogram(tmp,w1,ov,fftsize,250); % uses FFT
        tt2(k,:,:)=y;
    end
    tim=t;
    %compute power and plv
    %t=squeeze(abs(tt).^2);
    %t2=squeeze(mean(((t-repmat(mean(t,3),[1 1 size(y,2)]))./repmat(mean(t,3),[1 1 size(y,2)])),1));
    %allPOa(ic,:)=mean(mean(t(:,:,5:end-5),3),1);
    %allPO(ic,:,:)=squeeze(t2);
    allPL(ic,ic2,:,:)=squeeze(abs(mean(exp(i*(angle(tt)-angle(tt2))),1)));
    end
end
tim=tim-0.5;


 cfg=[];
 cfg.layout='4D248.lay';
 cfg.component=[1:20];
 ft_topoplotIC(cfg,comp);
 
%look at allPOa to get power spectra and find alpha components
alph=[1:3 8 12:15 19 20 22 25 28];

for k=alph,
    figure(1)
    pcolor(tim2,[1:40],squeeze(allPLdft(k,:,:)));shading interp;colorbar
    figure(2)
    pcolor(tim,f(1:80),squeeze(allPL(k,1:80,:)));shading interp;colorbar
    title(k)
    pause
end


for k=1:30,
    for k2=k+1:30,
    figure(1)
    pcolor(tim,f,squeeze(allPL(k,k2,:,:)));shading interp;colorbar
    title([k k2])
    pause
    end
end

2 14; 2 17; 2 20; 2 21; 3 26; 3 30; 4 22; 5 18;5 20; 6 8;6 18; 6 19; 8 16;