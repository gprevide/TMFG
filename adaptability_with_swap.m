addpath('C:\Users\Guido\Documents\Bauhaus\PhD\matlab_bgl', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\pmfg', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\BCT', ...
        'C:\Users\Guido\Documents\Bauhaus\PhD\RANDRAW');

%function adaptability_with_swap

clear 
close all
load('342USstocks.mat')
% there is one crazy price! YRC WORLDWIDE INC i=296 (but seems correct, it is an effect of the adjusting...)
figure,plot(dates,Prices(:,[1:295,297:end])),datetick 
Pr = Prices(:,[1:295,297:end]);
wi = 1000;
r = corrcoef(diff(log(Pr(1:wi,:))));
[P0, cliques, triangles0, tc0  ] = TMFGT1(r.^2);
P1=P0;P2=P0;P3=P0;
triangles1=triangles0;triangles2=triangles0;triangles3=triangles0;
tc1=tc0;tc2=tc0;tc3=tc0;
k = 0;

t=wi;
k=1;
%% from here can restart any time
t0=t;
k=k-1;

fileid = fopen('.\Paper_Tests_Real.csv', 'a');
fprintf(fileid, 'TMFG, TMFGT1, TMFG+T1, TMFG+T2, TMFGT2_K4, ORG+T2, ORG+T1, ORG+T2+T1, PMFG\n');
for t=t0:20:size(Pr,1)
    k = k+1;
    fprintf('timestep %d\n', t);
    t_d(k) = dates(t);
    r = corrcoef(diff(log(Pr((t-wi+1):t, :))));
    R = r.^2; R=(R+R')/2;
    P0(P0~=0) = R(P0~=0);
    [PT, cliquesT, trianglesT, tcT  ] = TMFG(R); % re-buid TMFG T2 all times
    [PT1, cliquesT1, trianglesT1, tcT1  ] = TMFGT1(R); % re-buid TMFG T1 all times
    [PTa, trianglesTa, tcTa] = T1(PT, trianglesT, tcT, R); PTa(PTa~=0)=R(PTa~=0);% add T1s
    [PTb, trianglesTb] = T2_Swap(PTa, R, trianglesTa); PTb(PTb~=0)=R(PTb~=0);% add Swaps
    [PTc, ~, ~, ~ ] = TMFGT2_K4(R); %last version of code where swaps are done while building the graph
    [P1, triangles1] = T2_Swap(P1, R, triangles1); P1(P1~=0)=R(P1~=0); % adding up swaps from original
    [P2, triangles2, tc2] = T1(P2, triangles2, tc2, R); P2(P2~=0)=R(P2~=0);% adding up T1 from original
    [P3, triangles3] = T2_Swap(P3, R, triangles3); P3(P3~=0)=R(P3~=0);
    [P3, triangles3, tc3] = T1(P3, triangles3, tc3, R); P3(P3~=0)=R(P3~=0); % adding up T1& Swaps from original
    PTP = pmfg(R); % re-buid PMFG all times
    sum_PT(k) = sum(PT(:));
    sum_PT1(k) = sum(PT1(:));
    sum_PTa(k) = sum(PTa(:));
    sum_PTb(k) = sum(PTb(:));    
    sum_PTc(k) = sum(PTc(:));
    sum_P0(k) = sum(P0(:));
    sum_P1(k) = sum(P1(:));
    sum_P2(k) = sum(P2(:));
    sum_P3(k) = sum(P3(:));
    sum_PTP(k) = sum(PTP(:));
    fprintf(fileid, '%f,%f,%f,%f,%f,%f,%f,%f,%f\n', sum_PT(k)+0.0, sum_PT1(k)+0.0, sum_PTa(k)+0.0, sum_PTb(k)+0.0, sum_PTc(k)+0.0, ...
        sum_P0(k)+0.0, sum_P1(k)+0.0, sum_P2(k)+0.0, sum_P3(k)+0.0, sum_PTP(k)+0.0);
        

    figure(3)
    clf
    plot(t_d(1:length(sum_PT)),sum_PT, 'k'); %TMFG T2
    hold on
    plot(t_d(1:length(sum_PT)),sum_PT1, '--k'); %TMFG T1
    plot(t_d(1:length(sum_PTa)),sum_PTa, ':k'); %TMFG+T1
    plot(t_d(1:length(sum_PTb)),sum_PTb, '-.k'); %TMFG+T1+Sw
    plot(t_d(1:length(sum_PTc)),sum_PTc, '-xk'); %'building with swaps
    plot(t_d(1:length(sum_P0)),sum_P0, '-y'); % Original
    plot(t_d(1:length(sum_P1)),sum_P1, 'r'); %swap only
    plot(t_d(1:length(sum_P2)),sum_P2, 'b'); %T1 only
    plot(t_d(1:length(sum_P3)),sum_P3, 'g'); %T1 only
    plot(t_d(1:length(sum_PTP)),sum_PTP, 'm'); %PMFG
    datetick
    legend({'TMFG T2','TMFG T2T1','TMFG+T1','TMFG+T1+Sw','building with swaps','original','swap only','T1 only','swap+T1','PMFG'},'location','northeast')
    drawnow
end

fclose(fileid);
%%
figure(2)
clf
%plot(t_d([1,end]),[0 0],'-k'); % zero
%hold on
plot(t_d(1:length(sum_PT1)),sum_PT1./sum_PT-1,'--k'); % TMFGT1
hold on
plot(t_d(1:length(sum_PTa)),sum_PTa./sum_PT-1,':k'); % TMFG+T1
plot(t_d(1:length(sum_PTb)),sum_PTb./sum_PT-1,'-.k'); %TMFG+Swap% 
plot(t_d(1:length(sum_PTc)),sum_PTc./sum_PT-1,'-xk'); %TMFG+v2% 
plot(t_d(1:length(sum_P1)),sum_P1./sum_PT-1,'r'); % swap only
plot(t_d(1:length(sum_P2)),sum_P2./sum_PT-1,'b'); % T1 only
plot(t_d(1:length(sum_P3)),sum_P3./sum_PT-1,'g'); % swap+T1
plot(t_d(1:length(sum_PT)),sum_PTP./sum_PT-1,'m'); % PMFG
plot(t_d(1:length(sum_P0)),sum_P0./sum_PT-1,'-y'); % original
datetick
drawnow
m=axis;
axis([m(1) m(2) -.3 0.02])
legend({'TMFGT1','TMFG+T1','TMFG+T1+Sw','TMFG-V2','swap only','T1 only','swap+T1','PMFG','original'},'location','southwest')
xlabel('time','fontsize',18,'interpreter','latex')
ylabel('relative score \% with respect to TMFG','fontsize',18,'interpreter','latex')
set(gca,'fontsize',12)
orient landscape
print('relative_scores.pdf','-dpdf')

%%
figure(3)
clf
plot(t_d(1:length(sum_PT)),sum_PT, 'k');
hold on
plot(t_d(1:length(sum_PT)),sum_PT1, '--k');
plot(t_d(1:length(sum_PT)),sum_PTa, ':k');
plot(t_d(1:length(sum_PT)),sum_PTb, '-.k');
plot(t_d(1:length(sum_PT)),sum_PTc, '-xk');
plot(t_d(1:length(sum_P1)),sum_P1, 'r');
plot(t_d(1:length(sum_P2)),sum_P2, 'b');
plot(t_d(1:length(sum_P3)),sum_P3, 'g');
plot(t_d(1:length(sum_PTP)),sum_PTP, 'm');
plot(t_d(1:length(sum_P0)),sum_P0, '-y');
datetick
m=axis;
axis([m(1) m(2) 70 260])
legend({'TMFG','TMFGT1','TMFG+T1','TMFG+T1+Sw','TMFG-V2','swap only','T1 only','swap+T1','PMFG','original'},'location','northwest')
xlabel('time','fontsize',18,'interpreter','latex')
ylabel('scores','fontsize',18,'interpreter','latex')
drawnow
set(gca,'fontsize',12)
orient landscape
print('scores.pdf','-dpdf')

%%
figure(4)
clf
%plot(t_d([1,end]),[0 0],'-k'); % zero
%hold on
[a,b]=hist(sum_PT./sum_PTP-1); % TMFGT1
plot(b,a,'-or')
hold on
[a,b]=hist(sum_PTa./sum_PTP-1); % TMFG+T1
plot(b,a,'-sb')
[a,b]=hist(sum_PTb./sum_PTP-1); %TMFG+Swap% 
plot(b,a,'-m^')
[a,b]=hist(sum_PTc./sum_PTP-1); %TMFG+v2% 
plot(b,a,'-kv')
legend({'TMFG T2','TMFG+T1','TMFG+T1+Sw','building with swaps'},'location','northeast')

% h = findobj(gca,'Type','patch'); 
% set(h(1),'FaceColor',[0.5 1 0.5],'facealpha',0.6)
% set(h(2),'FaceColor',[1 0 0],'facealpha',0.7)
% set(h(3),'FaceColor',[0 1 0],'facealpha',0.8)
% set(h(4),'FaceColor',[0 0 1],'facealpha',0.9)
% set(h(1),'EdgeColor',[0 0 0])


return
%%
stP=[];
sum_PTP=[];
%t0 = 2650;
t0=t;
r = corrcoef(diff(log(Pr((t0-wi+1):t0,:))));
% random start
[P0, cliques, triangles0, tc0  ] = TMFGT1(randn(size(r,1)).^2);
% start at point t0
% [P0, cliques, triangles0, tc0  ] = TMFGT1_mc(r.^2);
P3b = sparse(size(r,1));
P3b(P0~=0) = r(P0~=0).^2;
triangles3b = triangles0;
k=0;
for t=t0:20:size(Pr,1)
    k = k+1;
    fprintf('timestep %d\n', t);
    t_dP(k) = dates(t);
    r = corrcoef(diff(log(Pr((t-wi+1):t, :))));
    R = r.^2;
    PTP = doPMFG((R+R')/2); % re-buid TMFG T2 all times
    stP(k)=(sum((PTP(:)))-sum((PT(1:length(PTP)))))/sum((PTP(:))); % compute relative performancs difference
    sum_PTP(k) = sum(PTP(:));
    
    [P3b, triangles3b] = T2_Swap(P3b, R, triangles3b);
    P3b(P3b~=0)=R(P3b~=0);
    st3b(k)=(sum((P3b(:)))-sum((PT(1:length(PTP)))))/sum((P3b(:))); % compute relative performancs difference
    sum_P3b(k) = sum(P3b(:));

    
    figure(2)
    clf
    plot(t_d(1:length(st0)),st0,'--k'); 
    hold on
    plot(t_d(1:length(st1)),st1,'r'); 
    plot(t_dP(1:length(st3b)),st3b,':r'); 
    plot(t_d(1:length(st2)),st2,'b'); 
    plot(t_d(1:length(st3)),st3,'g'); 
    plot(t_dP(1:length(stP)),stP,'y'); 
    datetick
    drawnow

    figure(3)
    clf
    plot(t_d(1:length(sum_PT)),sum_PT, 'k');
    hold on
    plot(t_d(1:length(sum_P0)),sum_P0, '--k');
    plot(t_d(1:length(sum_P1)),sum_P1, 'r');
    plot(t_dP(1:length(sum_P3b)),sum_P3b, ':r');
    plot(t_d(1:length(sum_P2)),sum_P2, 'b');
    plot(t_d(1:length(sum_P3)),sum_P3, 'g');
    plot(t_dP(1:length(sum_PTP)),sum_PTP, 'y');
    datetick
    drawnow
end

