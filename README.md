# Economic-Dispatch-System
It is a system which will be used to reduce the fuel consumption of thermal power System.

#input data

import math

NG = 6              #number of generators=6
factor_delp=40.0    #factor_delp=cvp1(coefficient of variation of G=1)
factor_A=30.4        #A/sqr(delp)=2*a[1]
ITMAX=100           # maximum number of iterations
INMAX=100           #maximum number of iterations inner loop

a=[0.0,0.007,0.005,0.009,0.009,0.008,0.0075]
b=[0.0,7.0,10.0,8.5,11.0,10.5,12]
p=factor_A/(factor_delp*factor_delp)
print('p=%.9f'%p)
response=input('press <enter> key to continue:')

cvp=[0.0,0.012,0.04,0.012,0.012,0.012,0.012]
c=[0.0,240.0,200.0,220.0,200.0,220.0,120.0]

pgmin=[0.0,100.0,50.0,80.0,50.0,50.0,50.0]
pgmax=[0.0,500.0,200.0,300.0,150.0,200.0,120.0]

pd=700
alpha=0.005
beta=0.005
ep=0.001

B=[[0.056,0.0,0.0,0.0,0.0,0.0,0.0],
   [0.000,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001],
   [0.000,0.00001,0.00003,0.00001,0.00001,0.00001,0.00001],
   [0.000,0.00015,0.00001,0.00001,0.00001,0.00001,0.00001],
   [0.000,0.00001,0.00003,0.00001,0.00001,0.00001,0.00001],
   [0.000,0.00015,0.00001,0.00001,0.00001,0.00001,0.00001],
   [0.000,0.00001,0.00015,0.00001,0.00001,0.00001,0.00001]]

iter1=[]
gen1,gen2,la,gen3,gen4,gen5,gen6 = [],[],[],[],[],[],[]

for k in range (101): 
	iter1.append(0.0)
	la.append(0.0)
	gen1.append(0.0)
	gen2.append(0.0)
	gen3.append(0.0)
	gen4.append(0.0)
	gen5.append(0.0)
	gen6.append(0.0)
   
	
	
sum1=pd
sum2=0.0
for k in range(1,NG+1):
	sum1=sum1+(b[k]/(2*a[k]))
	sum2=sum2+(1/(2*a[k]))

lagr=sum1/sum2
print('initial value of lagrangian multiplier(lambda)=%.9f'%lagr)
response=input('press <enter> key to continue:')

pg=[[],[],[],[],[],[],[]]
for k in range (1001):
	pg[0].append(0.0)
	pg[1].append(0.0)
	pg[2].append(0.0)
	pg[3].append(0.0)
	pg[4].append(0.0)
	pg[5].append(0.0)
	pg[6].append(0.0)
                   
for k in range(1,NG+1):
	pg[k][0]=(lagr-b[k])/(2*a[k])
	print('initial value of p[%d]=%.7f'%(k,pg[k][0]))

pg[1][0]=0.0
pg[2][0]=0.0
pg[3][0]=0.0
pg[4][0]=0.0
pg[5][0]=0.0
pg[6][0]=0.0

u=0.05
ep_kt=56.145
OUTMAX=100
out=1

while(out<OUTMAX):          #outer loop
    it=1
    while(it<ITMAX):
        sum1=0.0
        in1=1
        P1,P2,P3,P4,P5,P6=0.0,0.0,0.0,0.0,0.0,0.0
        pg[1][it-1],pg[2][it-1],pg[3][it-1],pg[4][it-1],pg[5][it-1],pg[6][it-1]=0.0,0.0,0.0,0.0,0.0,0.0
        while(in1<INMAX):
            max=0.0
            for k in range(1,NG+1):
                sum1=0.0
                for m in range(1,NG+1):
                    sum1=sum1+2.0*B[k][m]*pg[m][it-1]
                pg[k][it]=(lagr*(1-sum1)-b[k])
                t2=(2*(a[k]+(a[k]*cvp[k]*cvp[k])+(p*cvp[k]*cvp[k])+(u*cvp[k]*cvp[k])+(lagr*B[k][k]*(1+cvp[k]*cvp[k]))))
                pg[k][it] = pg[k][it]/t2
                if(pg[k][it]<pgmin[k]):
                    pg[k][it]=pgmin[k]
                if(pg[k][it]>pgmax[k]):
                    pg[k][it]=pgmax[k]
            dif7=math.fabs(pg[1][it]-P1)
            if(dif7>max):
                max=dif7
            dif7=math.fabs(pg[2][it]-P2)
            if(dif7>max):
                max=dif7
            dif7=math.fabs(pg[3][it]-P3)
            if(dif7>max):
                max=dif7
            dif7=math.fabs(pg[4][it]-P4)
            if(dif7>max):
                max=dif7
            dif7=math.fabs(pg[5][it]-P5)
            if(dif7>max):
                max=dif7
            dif7=math.fabs(pg[6][it]-P6)
            if(dif7>max):
                max=dif7
                                  
            if(max<=0.001):
                print('it=%d,P1=%.7f,P2=%.7f lagr=%.9f'%(it,pg[1][it],pg[2][it],lagr))
                print('it=%d,P3=%.7f,P4=%.7f lagr=%.9f'%(it,pg[3][it],pg[4][it],lagr))
                print('it=%d,P5=%.7f,P6=%.7f lagr=%.9f'%(it,pg[5][it],pg[6][it],lagr))
                
                iter1[it]=it
                gen1[it]=pg[1][it]
                gen2[it]=pg[2][it]
                gen3[it]=pg[3][it]
                gen4[it]=pg[4][it]
                gen5[it]=pg[5][it]
                gen6[it]=pg[6][it]

                la[it]=lagr
                response=input('press <enter> key to continue:')
                break
            else:
                in1+=1
                P1=pg[1][it]
                P2=pg[2][it]
                P3=pg[3][it]
                P4=pg[4][it]
                P5=pg[5][it]
                P6=pg[6][it]
                pg[1][it-1]=pg[1][it]
                pg[2][it-1]=pg[2][it]
                pg[3][it-1]=pg[3][it]
                pg[4][it-1]=pg[4][it]
                pg[5][it-1]=pg[5][it]
                pg[6][it-1]=pg[6][it]
                continue
        sum1=B[0][0]
        for k in range (1,NG+1):
            sum1=sum1+(B[k][0]*pg[k][it]*pg[k][it])
        sum2=0.0
        for k in range (1,NG+1):
            for m in range (1,NG+1):
                sum2=sum2+(pg[k][it]*B[k][m]*pg[m][it])
        PL=sum1+sum2
        sump=0.0
        for k in range (1,NG+1):
            sump=sump+pg[k][it]
        deltap=(pd+PL-sump)
        deltap1=math.fabs(deltap)
        if(deltap1<=ep):
            print('success iteration = %d'%it)
            print('')
            optimal_total_cost=0.0
            for k in range (1,NG+1):
                print('IT=%d P=%d =%f lagr=%f deltap=%f ploss=%f'%(it,k,pg[k][it],lagr,deltap,PL))
                optimal_total_cost=optimal_total_cost+(a[k]*pg[k][it]*pg[k][it]+b[k]*pg[k][it]+c[k])

                if(k==1):
                    P1=pg[k][it]
                if(k==2):
                    P2=pg[k][it]
                if(k==3):
                    P3=pg[k][it]
                if(k==4):
                    P4=pg[k][it]
                if(k==5):
                    P5=pg[k][it]
                if(k==6):
                    P6=pg[k][it]
            print('optimal_total_cost=%f'%optimal_total_cost)
            break
        lagrnew=lagr+(alpha*deltap)
        lagr=lagrnew
        it+=1
    response=input('press <enter> key to continue:')
    diff1=cvp[1]*cvp[1]*P1*P1
    diff2=cvp[2]*cvp[2]*P2*P2
    diff3=cvp[3]*cvp[3]*P3*P3
    diff4=cvp[4]*cvp[4]*P4*P4
    diff5=cvp[5]*cvp[5]*P5*P5
    diff6=cvp[6]*cvp[6]*P6*P6
    diff=diff1+diff2+diff3+diff4+diff5+diff6
    print('it=%d,diff=%f,u=%f,lagr=%f'%(it,diff,u,lagr))
    print('P1=%f P2=%f P3=%f P4=%f P5=%f P6=%f'%(P1,P2,P3,P4,P5,P6))
    
    print('diff1=%f diff2=%f diff3=%f diff4=%f diff5=%f diff6=%f'%(diff1,diff2,diff3,diff4,diff5,diff6))
    
    response=input('press <enter> key to continue:')
    t1=math.fabs(diff-ep_kt)
    if(t1<=0.01):
        break
    else:
        u=u-(beta*(ep_kt-diff))
#outer loop closed
for k in range(1,11):
    print('it=%d P1=%f P2=%f lagr=%f'%(iter1[k],gen1[k],gen2[k],la[k]))
    print('it=%d P3=%f P4=%f lagr=%f'%(iter1[k],gen3[k],gen4[k],la[k]))
    print('it=%d P5=%f P6=%f lagr=%f'%(iter1[k],gen5[k],gen6[k],la[k]))
print('it=%d P1=%f P2=%f lagr=%f'%(it,pg[1][it],pg[2][it],lagr)) 
print('it=%d P3=%f P4=%f lagr=%f'%(it,pg[3][it],pg[4][it],lagr))
print('it=%d P5=%f P6=%f lagr=%f'%(it,pg[5][it],pg[6][it],lagr))
print('input data')
print('a1=%f b1=%f c1=%f cvp1=%f'%(a[1],b[1],c[1],cvp[1]))
print('a2=%f b2=%f c2=%f cvp2=%f'%(a[2],b[2],c[2],cvp[2]))
print('a3=%f b3=%f c3=%f cvp3=%f'%(a[3],b[3],c[3],cvp[3]))
print('a4=%f b4=%f c4=%f cvp4=%f'%(a[4],b[4],c[4],cvp[4]))
print('a5=%f b5=%f c5=%f cvp5=%f'%(a[5],b[5],c[5],cvp[5]))
print('a6=%f b6=%f c6=%f cvp6=%f'%(a[6],b[6],c[6],cvp[6]))

print('P1=%f '%P1)
print('P2=%f '%P2)
print('P3=%f '%P3)
print('P4=%f '%P4)
print('P5=%f '%P5)
print('P6=%f '%P6)


print('Risk (p) used in scheme=%f'%p)
print('p: A=%f  delp=%f (p= A/(delp*delp)'%(factor_A,factor_delp))
print('lagrangian multiplier=%f$/hr'%lagr)
print('optimal total cost=%f$/hr'%optimal_total_cost)
print('pd=%f(MW) Ploss=%f'%(pd,PL))
print('Variation=%f  u=%f'%(diff,u))
response=input('press <enter> key to continue:')
		
