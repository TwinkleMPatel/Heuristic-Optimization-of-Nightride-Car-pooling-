clear all
clc
clf
N=0;
  cities=20;% number of cities
  x=randi(100000,1,cities);
  x=sort(x);
  f=1;
  l1=2;
  pos=[];
  c1=[1 1];
  s=1;
  c2=[99 99];
  pos1=[];
  pos2=[];
  fz=1;
  fz1=1;
for time=1:1:100000
    if f<cities
       if time==x(f)
        f=f+1;
        k=100*rand(2,2);
        m=length(pos);
        pos(m+1,1)=k(1,1);
        pos(m+1,2)=k(1,2);
        pos(m+2,1)=k(2,1);
        pos(m+2,2)=k(2,2);
        N=N+1;
       end
        
    end
    
   if l1==length(pos)
       if isempty(pos1) && isempty(pos2)
          if dist1(c1,pos(m+1,:))<dist1(c2,pos(m+1,:))
         pos1(1,1)=pos(m+1,1);
         pos1(1,2)=pos(m+1,2);
         pos1(2,1)=pos(m+2,1);
         pos1(2,2)=pos(m+2,2);
         path1=pos1;
          else
              pos2(1,1)=pos(m+1,1);
              pos2(1,2)=pos(m+1,2);
              pos2(2,1)=pos(m+2,1);
              pos2(2,2)=pos(m+2,2);
              path2=pos2;
          end
      elseif isempty(pos1)
              pos1(1,1)=pos(m+1,1);
              pos1(1,2)=pos(m+1,2);
              pos1(2,1)=pos(m+2,1);
              pos1(2,2)=pos(m+2,2);
              path1=pos1;
      elseif isempty(pos2)
              pos2(1,1)=pos(m+1,1);
              pos2(1,2)=pos(m+1,2);
              pos2(2,1)=pos(m+2,1);
              pos2(2,2)=pos(m+2,2);
              path2=pos2;
      else % write else from here__________________________________
            lpos1=length(pos1);       
        for i=1:(lpos1/2)
            dpos1(i,1)=pos1((2*i),1);
            dpos1(i,2)=pos1((2*i),2);
        end
            dpos1(i+1,1)=pos(m+2,1);
            dpos1(i+1,2)=pos(m+2,2);

            lpos2=length(pos2);       
for i=1:(lpos2/2)
dpos2(i,1)=pos2((2*i),1);
dpos2(i,2)=pos2((2*i),2);
end
dpos2(i+1,1)=pos(m+2,1);
dpos2(i+1,2)=pos(m+2,2);

[bpos1,bdist1]=GA(length(dpos1)+1,c1,dpos1);% check lpos1 __________________
for u=1:length(dpos1)
    if bpos1(u,1)==pos(m+2,1)
        h1=u;
    end
end
[bpos2,bdist2]=GA(length(dpos2)+1,c2,dpos2);
for u=1:length(dpos2)
    if bpos2(u,1)==pos(m+2,1)
        h2=u;
    end
end

[odist1,path11]=mindist(c1,pos(m+1,:),bpos1,h1);
[odist2,path12]=mindist(c2,pos(m+1,:),bpos2,h2);

% change pos1 according to out1________________________________

if odist1<odist2
              pos1(lpos1+1,1)=pos(m+1,1);
              pos1(lpos1+1,2)=pos(m+1,2);
              pos1(lpos1+2,1)=pos(m+2,1);
              pos1(lpos1+2,2)=pos(m+2,2);
              path1=path11;
else
              pos2(lpos2+1,1)=pos(m+1,1);
              pos2(lpos2+1,2)=pos(m+1,2);
              pos2(lpos2+2,1)=pos(m+2,1);
              pos2(lpos2+2,2)=pos(m+2,2);
              path2=path12;
end
end
l1=l1+2;
   end
   
   if isempty(pos1)
   else
    if fz<length(path1)
     
  if sqrt(((c1(1,1)-path1(fz,1))^2)+((c1(2)-path1(fz,2))^2))>0.5
%       (abs(c1(1,1)-path1(i,1))>=0.5) || (abs(c1(2)-path1(i,2))>=0.5)
  theta=atan2((path1(fz,2)-c1(1,2)),(path1(fz,1)-c1(1,1)));
  c1=[(c1(1,1)+(cos(theta)*s)) (c1(1,2)+(sin(theta)*s))];

plot(path1(:,1),path1(:,2),'o');
hold on
% plot(c1(1),c1(2),'*r');
plot(c1(1),c1(2),'-s','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[0.8 0.9 1])
hold on
axis([0,100,0,100]);
drawnow
  end
  hold off
        if sqrt(((c1(1,1)-path1(fz,1))^2)+((c1(2)-path1(fz,2))^2))<=0.5
            if length(pos1)==2
                
            else
        fz=fz+1;
        for op=1:length(pos1)
            if mod(op,2)==0
            if path1(fz,1)==pos1(op,1)
                pos1(op,1)=0;
                pos1(op,2)=0;
                pos1(op-1,1)=0;
                pos1(op-1,2)=0;
            end
            end
        end
        for op=1:length(pos1)
            if mod(op,2)==0
               pos1=zero(pos1);
            end
        end
        disp('pos')
        disp(pos1)
        disp('path')
        disp(path1)
            end
        end
        
    end
    end
%    if isempty(pos2)
%        hold on
%    else
%   if fz1<length(path2)  
%   if sqrt(((c2(1,1)-path2(fz1,1))^2)+((c2(2)-path2(fz1,2))^2))>0.5
% %       (abs(c2(1,1)-path2(i1,1))>=0.5) || (abs(c2(2)-path2(i1,2))>=0.5)
%   theta=atan2((path2(fz1,2)-c2(1,2)),(path2(fz1,1)-c2(1,1)));
%   c2=[(c2(1,1)+(cos(theta)*s)) (c2(1,2)+(sin(theta)*s))];
% 
% plot(path2(:,1),path2(:,2),'*');
% hold on
% % plot(c2(1),c2(2),'*g');
% plot(c2(1),c2(2),'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .3])
% hold on
% axis([0,100,0,100]);
% drawnow
% 
%   end
%   hold off
%     if sqrt(((c2(1,1)-path2(fz1,1))^2)+((c2(2)-path2(fz1,2))^2))<=0.5
%     fz1=fz1+1;
%         for op=1:length(pos2)
%             if mod(op,2)==0
%             if path2(fz,1)==pos2(op,1)
%                 pos2(op,1)=0;
%                 pos2(op,2)=0;
%                 pos2(op-1,1)=0;
%                 pos2(op-1,2)=0;
%             end
%             end
%         end
%            for op=1:length(pos2)
%             if mod(op,2)==0
%                pos2=zero(pos2);
%             end
%            end
%       
%     end
%   end
%    end
end
%  plot(pos1(:,1),pos1(:,2),'*')
%    hold on 
%    plot(pos2(:,1),pos2(:,2),'o')
disp('end of sim')
% for i=1:1:100
%     plot(c(1,i),c(2,i),'*')
%     axis([-100,100,-100,100]);
%     drawnow
% end


%     for i=1:2:(2*N)
% plot(pos([i (i+1)],1),pos([i (i+1)],2))
% hold on
% drawnow
%     end
function [out,out2] = mindist(c,p,bpos,h)
 k=length(bpos);
  mat=[c(1) c(2);p(1) p(2)];
 for j=3:k+2
     mat(j,1)=bpos(j-2,1);
     mat(j,2)=bpos(j-2,2);
 end
for u=1:h
    dist(u)=0;
     for i=1:(length(mat)-1)
       dist(u)=sqrt((mat(i,1)-mat(i+1,1))^2+(mat(i,2)-mat(i+1,2))^2)+dist(u);
     end
     if dist(u)==min(dist)
     nmat=mat;
     end
     mat([u+1 u+2],1)=mat([u+2 u+1],1);
     mat([u+1 u+2],2)=mat([u+2 u+1],2);
end
 out= min(dist(u));
 for p=2:k+2
 out2(p-1,1)=nmat(p,1);
 out2(p-1,2)=nmat(p,2);
 end
 end
function [out,dist]= GA(N,c,pos1)
% clear all
% clc
% prompt='Enter value of N: ';
% N=input(prompt);
gen=50; % Total number of generations
x = ((N+3)/2)-mod((N+3)/2,1)-2; 
% d= N*(N+1)/2;
% x1=x+2;

posrt=[c(1);c(2)];
for j=2:N
    posrt(1,j)=pos1(j-1,1);
    posrt(2,j)=pos1(j-1,2);
end
% pos=100*rand(2,N); % assigning x and y values for N cities
th=2:N;
th=th(randperm(length(th)));
p1 =[1 th]; %_____________________________________________________________________change this line
th=th(randperm(length(th)));
p2 =[1 th];


for fd=1:1:gen  % for loop for number of generations
    
    % assigning values to positions for parent 1
    parent1=p1;
    parent2=p2;
for i=1:1:N
    for j=1:1:N
        if p1(i)==j
        parent1(1,i)=posrt(1,j);
        parent1(2,i)=posrt(2,j);
        end
    end
end
% parent1(1,N+1)=parent1(1,1);
% parent1(2,N+1)=parent1(2,1);


    % assigning values to positions for parent 2
for i=1:1:N
    for j=1:1:N
if p2(i)==j
    parent2(1,i)=posrt(1,j);
    parent2(2,i)=posrt(2,j);
end
    end
end
% parent2(1,N+1)=parent2(1,1);
% parent2(2,N+1)=parent2(2,1);

p3=p1;
p4=p2;

% swapping according to PMX
dummy=p1;
for i=x:1:x+2
p1(i)=p2(i);
end

for i=x:1:x+2
p2(i)=dummy(i);
end


% saving position of parent where we get repeated values after swapping for
% parent 1 and parent 2
if p1(x)==p2(x) && p1(x+1)==p2(x+1) && p1(x+2)==p2(x+2)
else
j=1;
k=0;
for i=1:1:N
if x~=i && x+1~=i && x+2~=i 
    if (p1(x)==p1(i) || p1(x+1)==p1(i) || p1(x+2)==p1(i))
    k(1,j)=i;
    if p1(x)==p1(i)
    k(2,j)=x;
    end
    if p1(x+1)==p1(i)
    k(2,j)=x+1;
    end
    if p1(x+2)==p1(i)
    k(2,j)=x+2;
    end
        j=j+1;
    end
end
end
k1=0;
j=1;
for i=1:1:N
if x~=i && x+1~=i && x+2~=i 
    if p2(x)==p2(i) || p2(x+1)==p2(i) || p2(x+2)==p2(i)
    k1(1,j)=i;
    if p2(x)==p2(i)
    k1(2,j)=x;
    end
    if p2(x+1)==p2(i)
    k1(2,j)=x+1;
    end
    if p2(x+2)==p2(i)
    k1(2,j)=x+2;
    end
    j=j+1;
    end
end
end
if k==0
    l=0;
else
l=length(k(1,:));
end
dummy3=p1;
if l==3
for g=1:1:3
p1(k(1,g))=p2(k(2,g));
end
for g1=1:1:3
    p2(k1(1,g1))=dummy3(k1(2,g1));
end
end
if l==2
        p1(k(1,1))=p2(k1(1,2));
        p1(k(1,2))=p2(k1(1,1));
        p2(k1(1,2))=dummy3(k(1,1));
        p2(k1(1,1))=dummy3(k(1,2));
end
    if l==1
        p1(k(1,1))=p2(k1(1,1));
        p2(k1(1,1))=dummy3(k(1,1));
    end
end


% mutation with probabity 
prob=0.8; % probability of mutation
if rand(1)<prob
% mut=randperm(N);
p1(1,[N (N-1)])=p1(1,[(N-1) (N)]);
end
if rand(1)<prob
% mut=randperm(N);
p1(1,[N (N-1)])=p1(1,[(N-1) (N)]);
end


% assigning position for child 1 and child 2
child1=p1;
child2=p2;
for i=1:1:N
    for j=1:1:N
if p1(i)==j
    child1(1,i)=posrt(1,j);
    child1(2,i)=posrt(2,j);
end
    end
end
% child1(1,N+1)=child1(1,1);
% child1(2,N+1)=child1(2,1);
for i=1:1:N
    for j=1:1:N
if p2(i)==j
    child2(1,i)=posrt(1,j);
    child2(2,i)=posrt(2,j);
end
    end
end
% child2(1,N+1)=child2(1,1);
% child2(2,N+1)=child2(2,1);



% calculating distance for all parents and children
pd1=distance(parent1,N);
pd2=distance(parent2,N);
cd1=distance(child1,N);
cd2=distance(child2,N);


% selecting 4 best values

bm=[pd1 pd2 cd1 cd2];... pd11 pd12 cd3 cd4];
temp = [p3;p4;p1;p2];...;p13;p14;p11;p12];
% disp('minimum distance among all parents and child:')
% min(bm)
bm1=sort(bm);
np=[];
dist=min(bm);
for i = 1:4
    for j = 1:4
        if bm1(i)==bm(j)
           np(i,:)=temp(j,:);
        end
    end    
end

p1=np(1,:);
% for plotting graph
for i=1:1:N
    for j=1:1:N
        if p1(i)==j
        bparent(1,i)=posrt(1,j);
        bparent(2,i)=posrt(2,j);
        end
    end
end
% bparent(1,N+1)=bparent(1,1);
% bparent(2,N+1)=bparent(2,1);
% plot(bparent(1,:),bparent(2,:),'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6])
% pause(0.05);
p2=np(2,:);

for j=2:N
    out(j-1,1)=bparent(1,j);
    out(j-1,2)=bparent(2,j);

end
end
end

function D=distance(L,N)
D=0;
   for j=1:1:(N-1)
   D=sqrt((L(1,j)-L(1,j+1))^2+(L(2,j)-L(2,j+1))^2)+D;
   end
end

function d=dist1(T,U)
d=sqrt((T(1,1)-U(1,1))^2 +(T(1,2)-U(1,2))^2);
end
function out= zero(a)
k=length(a);
if k==2
    out=[];
else
for i=1:k-1
    if i<=k
while a(i,1)==0 
    for q=i:(length(a)-1)
        a([q q+1],1)=a([q+1 q],1);
        a([q q+1],2)=a([q+1 q],2);
    end
    for q1=1:(length(a)-1)
        a1(q1,1)=a(q1,1);
        a1(q1,2)=a(q1,2);
    end
    a=a1;
    k=length(a);
    a1=[];
     if min(min(a))~=0
     break
     end
end
    end
end
out=a;
end
end
 