clear
clc
num = [400:527];
prefix_name = {'/media/user/9EAEE48CAEE45DF1/UQ/bed/'};
xx = 200:10:540;
zz = 1:65;
ttt = 1e-4;
lll = 0.038/10;
vvv = lll/ttt; % hello chenyilin how are you ?
% donnt be afraid, I am not ghost!!! I just study in  Jin's office for a
% few days.
% why are you often use Mathlab? I think C++ would be more efficient.
%who are you?
%C++ is most used for simulatio this post process
%Are you Zhongtian shixiong?
%do not  turn off the computer, still running the simulation
% yes, ok, i won't turn off it
for i = 1:numel(num)
    name = strcat(prefix_name,'perbed__',num2str(num(i),'%04d'),'.h5');
	[data, domain] = getData(char(name),char('/Density_0'),2,false);
	Nx = domain.Nx;
	Ny = domain.Ny;
	Nz = domain.Nz;
    rho = reshape(data(:,:,:),[Nx,Ny,Nz]);
    p = rho./3;
    px = sum(sum(p(xx,2:Ny-1,2:Nz-1),3),2);
    subplot(2, 1, 1)
    plot(xx,px,'-*')
    title(strcat(num2str(num(i)),' s'))
    kk = polyfit(xx',px,1);
    hold on
    pxf = kk(1).*xx+kk(2);
    plot(xx,pxf,'r')
    drawnow	
    hold off
    g(i,1) = kk(1);
    subplot(2, 1, 2)
    plot(i+num(1),kk(1),'-s')
    hold on
end %


