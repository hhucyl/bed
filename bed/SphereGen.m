clear
clc
prefix_name = {'td_more2_'};
num = 9;
i =1;

name = strcat(prefix_name,num2str(num(i),'%04d'),'.h5');
post = getData(char(name),char('/Position'),1,false);
rt = getData(char(name),char('/Radius'),1,false);
kkk = find(rt ~= 0);
pos(:,1) = post(3.*kkk-2);
pos(:,2) = post(3.*kkk-1);
pos(:,3) = post(3.*kkk);
nn = 1;
xlim = [0 140];
ylim = [0 50];
zlim = [0 60];
for ip = 1:numel(pos(:,1))
    flag = pos(ip,1)>=xlim(1) & pos(ip,1)<=xlim(2) & pos(ip,2)>=ylim(1)...
        &pos(ip,2)<=ylim(2) & pos(ip,3)>=zlim(1) & pos(ip,3)<=zlim(2);
    if(flag)
        ppos(nn,:) = pos(ip,:) - [xlim(1),ylim(1),0];
        nn = nn+1;
    end
end
for i=1:nn-1;
    sphereplot(ppos(i,1),ppos(i,2),ppos(i,3),5)
    drawnow
end
fid = fopen('sphere.txt','w');
fprintf(fid,'%d\n',nn-1);
fprintf(fid,'%f %f %f\n',ppos');
fclose(fid)

