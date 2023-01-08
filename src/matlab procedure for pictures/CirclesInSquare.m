clear
clc
filename = '.\153_Sol.txt';
[x,y] = textread(filename,'%n%n') ;
Ncircles = x(1) ;
R = 0.5*y(1) ; 
Circle  = zeros(Ncircles,2) ;
contact = zeros(Ncircles,1);
for i=1:Ncircles
    Circle(i,1) = x(i+1);
    Circle(i,2) = y(i+1); 
end 
figure;
hold on;
line([-1.0*R,R],[R,R], 'linewidth',1.3, 'linestyle','-','color','k');
line([R,R],[R,-1.0*R], 'linewidth',1.3, 'linestyle','-','color','k');
line([R,-1.0*R],[-1.0*R,-1.0*R], 'linewidth', 1.3,'linestyle','-','color','k');
line([-1.0*R,-1.0*R],[-1.0*R,R], 'linewidth', 1.3,'linestyle','-','color','k');

axis equal
axis off
epslon = 1.0e-6; 
for i = 1:Ncircles
    r1 = abs(Circle(i,1));
    if r1 + 1.0 > R - epslon
       contact(i) =  contact(i)+1;  
    end
end
for i = 1:Ncircles
    r1 = abs(Circle(i,2)); 
    if r1 + 1.0 > R - epslon
       contact(i) =  contact(i)+1;  
    end
end

for i = 1:Ncircles-1
    for j = i+1:Ncircles
        if sqrt((Circle(i,1)-Circle(j,1))*(Circle(i,1)-Circle(j,1)) + (Circle(i,2)-Circle(j,2))*(Circle(i,2)-Circle(j,2)))  < 2.0 + epslon
          contact(i) =  contact(i)+1;
          contact(j) =  contact(j)+1;
        end 
    end
end

t= 0:0.001:2.0*pi;
for i =1:Ncircles
    c1= 1.0*cos(t)+Circle(i,1);
    c2= 1.0*sin(t)+Circle(i,2); 
    plot(c1,c2,'-k'); 
    switch contact(i)
    case 0
         fill(c1,c2,[1 0 0]); %[1 0 0] is RGB of color
    case 1
          fill(c1,c2,[255/255.0 105/255.0 180/255.0]);   
    case 2
          fill(c1,c2,[0 191/255.0 255/255.0]); 
    case  3
          fill(c1,c2,[0.667 0.667 1]);   
    case  4
          fill(c1,c2,'g');  
    case  5
          fill(c1,c2,[238/255.0 221/255.0 130/255.0]);   
    otherwise 
          fill(c1,c2,[255/255 165/255 0]);          
    end 
end 

for i = 1:Ncircles-1
    for j = i+1:Ncircles
        if sqrt((Circle(i,1)-Circle(j,1))*(Circle(i,1)-Circle(j,1)) + (Circle(i,2)-Circle(j,2))*(Circle(i,2)-Circle(j,2)))  < 2.0 + epslon
          line([Circle(i,1),Circle(i,1)+0.3*(Circle(j,1)-Circle(i,1))],[Circle(i,2),Circle(i,2)+0.3*(Circle(j,2)-Circle(i,2))],'linestyle','-','color','k');
          line([Circle(j,1),Circle(j,1)+0.3*(Circle(i,1)-Circle(j,1))],[Circle(j,2),Circle(j,2)+0.3*(Circle(i,2)-Circle(j,2))],'linestyle','-','color','k');
        end
    end
end

for i = 1:Ncircles
    r1 = abs(Circle(i,1));
    if r1 > R - 1.0 - epslon
        if Circle(i,1) > 0
           line([Circle(i,1),Circle(i,1)+0.5],[Circle(i,2), Circle(i,2)],'linestyle','-','color','k'); 
        else
           line([Circle(i,1),Circle(i,1)-0.5],[Circle(i,2), Circle(i,2)],'linestyle','-','color','k');   
        end
    end
end
for i = 1:Ncircles
    r1 = abs(Circle(i,2));
    if r1 > R - 1.0 - epslon
        if Circle(i,2) > 0
           line([Circle(i,1),Circle(i,1)],[Circle(i,2), Circle(i,2)+0.5],'linestyle','-','color','k'); 
        else
           line([Circle(i,1),Circle(i,1)],[Circle(i,2), Circle(i,2)-0.5],'linestyle','-','color','k');   
        end
    end
end


