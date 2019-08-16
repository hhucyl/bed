clear
clc
num = [496];
prefix_name = {'/media/user/9EAEE48CAEE45DF1/UQ/bed/'};
xx = 200:10:400;
zz = 1:65;
ttt = 1e-4;
lll = 0.038/10;
vvv = lll/ttt;
exdata = [ 1.20132013201320e-001    4.81189851268589e-002
 1.17491749174918e-001    3.93700787401574e-002
 1.14888155482215e-001    2.18722659667536e-002
 1.10231023102310e-001    0.00000000000000e+000
 9.47194719471947e-002   -5.90551181102361e-002
 7.45874587458746e-002   -1.11548556430447e-001
 9.64429776310964e-003   -9.90813648293964e-001
 1.53648698203154e-002   -1.99037620297463e+000
 1.55115511551155e-002   -2.14785651793526e+000
 1.65016501650165e-002   -2.99212598425197e+000
 1.71617161716172e-002   -3.14523184601925e+000
 2.11221122112211e-002   -3.97637795275591e+000
 2.47158049138247e-002   -4.98468941382327e+000];

for i = 1:numel(num)
    name = strcat(prefix_name,'perbed__',num2str(num(i),'%04d'),'.h5');
	[data, domain] = getData(char(name),char('/Velocity_0'),2,true);
	Nx = domain.Nx;
	Ny = domain.Ny;
	Nz = domain.Nz;
    vx = reshape(data(:,:,:,1),[Nx,Ny,Nz]);
    vy = reshape(data(:,:,:,2),[Nx,Ny,Nz]);
    vz = reshape(data(:,:,:,3),[Nx,Ny,Nz]);
    vvx = sum(sum(vx(xx,2:Ny-1,:),2),1)/((Ny-2)*numel(xx));
    vvx = vvx(:);
    plot(exdata(:,1)./vvv,(exdata(:,2)./(lll/0.038))./10,'ro')
    hold on
    plot(vvx(zz),zz./10-6.2,'*')
    title(strcat(num2str(num(i)),' s'))
    ylabel('\itz^*')
    xlabel('u m/s')
    legend('Experiment Re = 39000','Stimulation','location','Southeast')
    hold off
    drawnow	 
end %