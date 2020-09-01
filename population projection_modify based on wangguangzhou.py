import math

global axn_m
global axn_f
axn_m = list()
axn_f = list()

def male_ax(dr0,pk):
    age_mqx = 1- math.exp(-dr0)
    if (age_mqx > 0.100):
        axn_m[0]=0.33
        for i in range(1,5):
            axn_m[i] = 1.352/4
        for i in range(5,10):
            axn_m[i]=2.25/5
        for i in range(10,71):
            axn_m[i]=2.6/5
        for i in range(71,pk):
            axn_m[i]=2.5/5
    else:
        axn_m[0]=0.0425 + 2.875*age_mqx
        for i in range(1,5):
            axn_m[i]=(1.63-3.013*age_mqx)/4
        for i in range(5,10):
            axn_m[i]=2.25/5
        for i in range(10,71):
            axn_m[i]=2.6/5
        for i in range(71,pk):
            axn_m[i]=2.5/5
    return axn_m



def female_ax(dr0,pk):
    age_fqx = 1- math.exp(-dr0)
    if (age_fqx > 0.100):
        axn_f[0]=0.35
        for i in range(1,5):
            axn_f[i] = 1.364/4
        for i in range(5,10):
            axn_f[i]=2.25/5
        for i in range(10,71):
            axn_f[i]=2.6/5
        for i in range(71,pk):
            axn_f[i]=2.5/5
    else:
        axn_f[0]=0.0425 + 2.875*age_fqx
        for i in range(1,5):
            axn_f[i]=(1.524-1.627*age_fqx)/4
        for i in range(5,10):
            axn_f[i]=2.25/5
        for i in range(10,71):
            axn_f[i]=2.6/5
        for i in range(71,pk):
            axn_f[i]=2.5/5
    return axn_f



def ex_m(dr,pk):
    mqx=list()
    mdx=list()
    mlx=list()
    mCLx=list()
    mTx=list()
    mex=list()
    for i in range(pk):
        mqx.append(0)
        mdx.append(0)
        mlx.append(0)
        mCLx.append(0)
        mTx.append(0)
        mex.append(0)
        
    mlx.append(0)
    for i in range(pk-1):
        mqx[i]=1-math.exp(-dr[i])
        
    mqx[pk-1]=1
    mlx[0]=1
    for i in range(pk):
        mdx[i]=mlx[i]*mqx[i]
        mlx[i+1]=mlx[i]-mdx[i]
    for i in range(0,pk-1):
        mCLx[i]=mlx[i+1]-axn_m[i]*(mlx[i]-mlx[i+1])
    mCLx[pk-1]=mlx[pk-1]/dr[pk-1]
    for i in range(pk):
        temp=0
        for j in range(i,pk):
            temp=temp+mCLx[j]
            mTx[i]=temp
        mex[i]=mTx[i]/mlx[i]
    return mlx,mCLx,mex



def ex_f(dr,pk):
    fqx=list()
    fdx=list()
    flx=list()
    fCLx=list()
    fTx=list()
    fex=list()
    for i in range(pk):
        fqx.append(0)
        fdx.append(0)
        flx.append(0)
        fCLx.append(0)
        fTx.append(0)
        fex.append(0)
        
    flx.append(0)
    
    for i in range(pk-1):
        fqx[i]=1-math.exp(-dr[i])
        
    fqx[pk-1]=1
    flx[0]=1
    for i in range(pk):
        fdx[i]=flx[i]*fqx[i]
        flx[i+1]=flx[i]-fdx[i]
    for i in range(0,pk-1):
        fCLx[i]=flx[i+1]-axn_f[i]*(flx[i]-flx[i+1])
    fCLx[pk-1]=flx[pk-1]/dr[pk-1]
    for i in range(pk):
        temp=0
        for j in range(i,pk):
            temp=temp+fCLx[j]
            fTx[i]=temp
        fex[i]=fTx[i]/flx[i]
    return flx,fCLx,fex



def lx_male(mex0,dr_e,mex_p,mlx_p,pk):
    epslon=0.000001
    mlx=list()
    mCLx=list()
    mTx=list()
    mex=list()
    for i in range(pk):
        mlx.append(0)
        mCLx.append(0)
        mTx.append(0)
        mex.append(0)
    
    mlx.append(0)
    
    arfa_min=-100.00
    arfa_max=100.00
    current_e0=0.0
    male_arfa=0.0
    keys=100.00
    
    while(keys>=epslon):
        male_arfa=(arfa_max+arfa_min)/2
        current_e0=0.0
        mlx_p[0]=1
        for i in range(1,pk):
            mlx[i]=mlx_p[i]/(mlx_p[i]+(1-mlx_p[i])*math.exp(2*male_arfa))
        mlx[0]=1
        for i in range(0,pk-1):
            mCLx[i]=mlx[i+1]+axn_m[i]*(mlx[i]-mlx[i+1])
        # mCLx[pk-1]=mlx[pk-1]/dr_e
        # mCLx[pk-1]=mlx[pk-1]+(mlx[pk-1]-mlx[pk])/2
        for i in range(pk):
            temp = 0
            for j in range(i,pk):
                temp=temp+mCLx[j]
                mTx[i]=temp
            mex[i]=mTx[i]/mlx[i]
        current_e0=mex[0]
        if current_e0<mex_p:
            arfa_max=male_arfa
            keys= abs(mex_p-current_e0)
        else:
            arfa_min=male_arfa
            keys=abs(mex_p-current_e0)
        return mCLx



def lx_female(fex0,dr_e,fex_p,flx_p,pk):
    epslon=0.000001
    flx=list()
    fCLx=list()
    fTx=list()
    fex=list()
    for i in range(pk):
        flx.append(0)
        fCLx.append(0)
        fTx.append(0)
        fex.append(0)
    
    flx.append(0)
    
    arfa_min=-100.00
    arfa_max=100.00
    current_e0=0.0
    female_arfa=0.0
    keys=100.00
    
    while(keys>=epslon):
        female_arfa=(arfa_max+arfa_min)/2
        current_e0=0.0
        flx_p[0]=1
        for i in range(1,pk):
            flx[i]=flx_p[i]/(flx_p[i]+(1-flx_p[i])*math.exp(2*female_arfa))
        flx[0]=1
        for i in range(0,pk-1):
            fCLx[i]=flx[i+1]+axn_f[i]*(flx[i]-flx[i+1])
        # fCLx[pk-1]=flx[pk-1]/dr_e
        fCLx[pk-1]=flx[pk-1]+(flx[pk-1]-flx[pk])/2
        for i in range(pk):
            temp = 0
            for j in range(i,pk):
                temp=temp+fCLx[j]
                fTx[i]=temp
            fex[i]=fTx[i]/flx[i]
        current_e0=fex[0]
        if current_e0<fex_p:
            arfa_max=female_arfa
            keys= abs(fex_p-current_e0)
        else:
            arfa_min=female_arfa
            keys=abs(fex_p-current_e0)
        return fCLx




infile=open('/Users/apple/desktop/usualdata/chancheng_changzhu_population_data.csv','r')
inparm=open('/Users/apple/desktop/usualdata/chancheng_changzhu_population_project.csv','r')
outfile=open('/Users/apple/desktop/usualdata/chancheng_changzhu_population_output.csv','w')
outfile2=open('/Users/apple/desktop/usualdata/chancheng_changzhu_population_output2.csv','w')

age=list()
param=list()

mdr=list()
fdr=list()

male=list()
female=list()
male_lx=list()
female_lx=list()
male_blx=list()
female_blx=list()
asfr_all=list()

mex_p=list()
fex_p=list()

Tfr_p=list()
SRB=list()
Tbirth=list()

Tpop=list()
T0_14=list()
T15_59=list()
T15_64=list()
T60plus=list()
T65plus=list()

T0_4=list()
T5_9=list()
T10_14=list()
T15_19=list()
T20_24=list()
T25_29=list()
T30_34=list()
T35_39=list()
T40_44=list()
T45_49=list()
T50_54=list()
T55_59=list()
T60_64=list()
T65_69=list()
T70_74=list()
T75_79=list()
T80_84=list()
T85_89=list()
T90_94=list()
T95_99=list()

temp_male=list()
temp_female=list()
temp_male_blx=list()
temp_female_blx=list()
temp_asfr=list()



for row in infile:
    aa=row.split(',')
    pj=len(aa)
    age.append(aa)

pl=len(age)
pk=pl-1

for i in range(pk):
    axn_m.append(0)
    axn_f.append(0)

for i in range(pk):
    temp_male.append(float(age[i+1][1]))
    temp_female.append(float(age[i+1][2]))
    mdr.append(float(age[i+1][3])/1000)
    fdr.append(float(age[i+1][4])/1000)
    temp_asfr.append(float(age[i+1][5])/1000)
    temp_male_blx.append(0)
    temp_female_blx.append(0)
    
for row in inparm:
    bb=row.split(',')
    py=len(bb)
    param.append(bb)
    
pt=len(param)
project_years=pt-1



for i in range(project_years):
    male.append([])
    female.append([])
    male_lx.append([])
    female_lx.append([])
    male_blx.append([])
    female_blx.append([])
    asfr_all.append([])
    
    mex_p.append([])
    fex_p.append([])
    Tfr_p.append([])
    SRB.append([])
    Tbirth.append(0.0)
    Tpop.append(0.0)
    T0_14.append(0.0)
    T15_59.append(0.0)
    T15_64.append(0.0)  
    T60plus.append(0.0)
    T65plus.append(0.0)

    T0_4.append(0.0)
    T5_9.append(0.0)
    T10_14.append(0.0)
    T15_19.append(0.0)
    T20_24.append(0.0)
    T25_29.append(0.0)
    T30_34.append(0.0)
    T35_39.append(0.0)
    T40_44.append(0.0)
    T45_49.append(0.0)
    T50_54.append(0.0)
    T55_59.append(0.0)
    T60_64.append(0.0)
    T65_69.append(0.0)
    T70_74.append(0.0)
    T75_79.append(0.0)
    T80_84.append(0.0)
    T85_89.append(0.0)
    T90_94.append(0.0)
    T95_99.append(0.0)
    
    
    for j in range(pk):
        male[i].append([])
        female[i].append([])
        male_lx[i].append([])
        female_lx[i].append([])
        male_blx[i].append([])
        female_blx[i].append([])
        asfr_all[i].append([])



for i in range(project_years):
    mex_p[i]=float(param[i+1][1])
    fex_p[i]=float(param[i+1][2])
    Tfr_p[i]=float(param[i+1][3])
    SRB[i]=float(param[i+1][4])


(axn_m) = male_ax(mdr[0], pk)
(mlx,mCLx,mex) = ex_m(mdr,pk)
(axn_f) = female_ax(fdr[0], pk)
(flx,fCLx,fex) = ex_f(fdr, pk)

temp_Tfr = 0
for j in range(pk):
    temp_Tfr = temp_Tfr+temp_asfr[j]

for j in range(pk):
    male[0][j] = temp_male[j]
    female[0][j] = temp_female[j]
    male_lx[0][j] = mlx[j]
    female_lx[0][j] = flx[j]
    male_blx[0][j] = mCLx[j]
    female_blx[0][j] = fCLx[j]
    asfr_all[0][j] = temp_asfr[j]/temp_Tfr

for i in range(1,project_years):
    for j in range(pk):
        asfr_all[i][j] = asfr_all[0][j]*Tfr_p[i]

for i in range(1,project_years):
    temp_male_blx = lx_male(mex_p[0], mdr[pk-1], mex_p[i], mlx, pk)
    for j in range(pk):
        male_blx[i][j]=temp_male_blx[j]

for i in range(1,project_years):
    temp_female_blx = lx_female(fex_p[0], fdr[pk-1], fex_p[i], flx, pk)
    for j in range(pk):
        female_blx[i][j]=temp_female_blx[j]



Tbirth[0]=male[0][0]+female[0][0]
for j in range(pk):
    Tpop[0]=Tpop[0]+male[0][j]+female[0][j]
for j in range(0,15):
    T0_14[0]=T0_14[0]+male[0][j]+female[0][j]
for j in range(15,60):
    T15_59[0]=T15_59[0]+male[0][j]+female[0][j]
for j in range(15,65):
    T15_64[0]=T15_64[0]+male[0][j]+female[0][j]
for j in range(60,pk):
    T60plus[0]=T60plus[0]+male[0][j]+female[0][j]
for j in range(65,pk):
    T65plus[0]=T65plus[0]+male[0][j]+female[0][j]  


for j in range(0,5):
    T0_4[0]=T0_4[0]+male[0][j]+female[0][j]
for j in range(5,10):
    T5_9[0]=T5_9[0]+male[0][j]+female[0][j]
for j in range(10,15):
    T10_14[0]=T10_14[0]+male[0][j]+female[0][j]

for j in range(15,20):
    T15_19[0]=T15_19[0]+male[0][j]+female[0][j]
for j in range(20,25):
    T20_24[0]=T20_24[0]+male[0][j]+female[0][j]
for j in range(25,30):
    T25_29[0]=T25_29[0]+male[0][j]+female[0][j]
for j in range(30,35):
    T30_34[0]=T30_34[0]+male[0][j]+female[0][j]
for j in range(35,40):
    T35_39[0]=T35_39[0]+male[0][j]+female[0][j]
for j in range(40,45):
    T40_44[0]=T40_44[0]+male[0][j]+female[0][j]
for j in range(45,50):
    T45_49[0]=T45_49[0]+male[0][j]+female[0][j]
for j in range(50,55):
    T50_54[0]=T50_54[0]+male[0][j]+female[0][j]
for j in range(55,60):
    T55_59[0]=T55_59[0]+male[0][j]+female[0][j]

for j in range(60,65):
    T60_64[0]=T60_64[0]+male[0][j]+female[0][j]
for j in range(65,70):
    T65_69[0]=T65_69[0]+male[0][j]+female[0][j]
for j in range(70,75):
    T70_74[0]=T70_74[0]+male[0][j]+female[0][j]
for j in range(75,80):
    T75_79[0]=T75_79[0]+male[0][j]+female[0][j]
for j in range(80,85):
    T80_84[0]=T80_84[0]+male[0][j]+female[0][j]
for j in range(85,90):
    T85_89[0]=T85_89[0]+male[0][j]+female[0][j]
for j in range(90,95):
    T90_94[0]=T90_94[0]+male[0][j]+female[0][j]
for j in range(95,100):
    T95_99[0]=T95_99[0]+male[0][j]+female[0][j]



for i in range(1,project_years):
    for j in range(pk-1):
        male[i][j+1]=male[i-1][j]*(male_blx[i][j+1]/male_blx[i][j])
        female[i][j+1]=female[i-1][j]*(female_blx[i][j+1]/female_blx[i][j])
    temp0 = 0.0
    for j in range(15,pk):
        temp0=temp0+0.5*(female[i-1][j]+
                         female[i-1][j-1]*female_blx[i][j]/female_blx[i][j-1])*asfr_all[i][j]
        Tbirth[0]=temp0
        male[i][0]=temp0*(SRB[i]/(100+SRB[i]))*male_blx[i][0]
        female[i][0]=temp0*(1-SRB[i]/(100+SRB[i]))*female_blx[i][0]
        
    for j in range(pk):
        Tpop[i]=Tpop[i]+male[i][j]+female[i][j]
    for j in range(0,1):
        Tbirth[i]=Tbirth[i]+male[i][j]+female[i][j]
    for j in range(0,15):
        T0_14[i]=T0_14[i]+male[i][j]+female[i][j]
    for j in range(15,60):
        T15_59[i]=T15_59[i]+male[i][j]+female[i][j]
    for j in range(15,65):
        T15_64[i]=T15_64[i]+male[i][j]+female[i][j]
    for j in range(60,pk):
        T60plus[i]=T60plus[i]+male[i][j]+female[i][j]
    for j in range(65,pk):
        T65plus[i]=T65plus[i]+male[i][j]+female[i][j]   
        
    for j in range(0,5):
        T0_4[i]=T0_4[i]+male[i][j]+female[i][j]
    for j in range(5,10):
        T5_9[i]=T5_9[i]+male[i][j]+female[i][j]
    for j in range(10,15):
        T10_14[i]=T10_14[i]+male[i][j]+female[i][j]
    for j in range(15,20):
        T15_19[i]=T15_19[i]+male[i][j]+female[i][j]
    for j in range(20,25):
        T20_24[i]=T20_24[i]+male[i][j]+female[i][j]
    for j in range(25,30):
        T25_29[i]=T25_29[i]+male[i][j]+female[i][j]
    for j in range(30,35):
        T30_34[i]=T30_34[i]+male[i][j]+female[i][j]
    for j in range(35,40):
        T35_39[i]=T35_39[i]+male[i][j]+female[i][j]
    for j in range(40,45):
        T40_44[i]=T40_44[i]+male[i][j]+female[i][j]
    for j in range(45,50):
        T45_49[i]=T45_49[i]+male[i][j]+female[i][j]
    for j in range(50,55):
        T50_54[i]=T50_54[i]+male[i][j]+female[i][j]
    for j in range(55,60):
        T55_59[i]=T55_59[i]+male[i][j]+female[i][j]
    for j in range(60,65):
        T60_64[i]=T60_64[i]+male[i][j]+female[i][j]
    for j in range(65,70):
        T65_69[i]=T65_69[i]+male[i][j]+female[i][j]
    for j in range(70,75):
        T70_74[i]=T70_74[i]+male[i][j]+female[i][j]
    for j in range(75,80):
        T75_79[i]=T75_79[i]+male[i][j]+female[i][j]
    for j in range(80,85):
        T80_84[i]=T80_84[i]+male[i][j]+female[i][j]
    for j in range(85,90):
        T85_89[i]=T85_89[i]+male[i][j]+female[i][j]
    for j in range(90,95):
        T90_94[i]=T90_94[i]+male[i][j]+female[i][j]
    for j in range(95,100):
        T95_99[i]=T95_99[i]+male[i][j]+female[i][j]


outfile.write('Year,Tpop,Tbirth,0-14,15-59,15-64,60+,65+,0-4,5-9,10-14,15-19,20-24,\
              25-29,30-34,35-39,40-44,45-49,50-54,55-59,60-64,65-69,70-74,75-79,\
              80-84,85-89,90-94,95-99')
outfile.write('\n')

for i in range(project_years):
    Lifetable = str(i)+','+str(Tpop[i])+','+str(Tbirth[i])+','+str(T0_14[i])+','\
            +str(T15_59[i])+','+str(T15_64[i])+','+str(T60plus[i])+','+str(T65plus[i])+','\
            +str(T0_4[i])+','+str(T5_9[i])+','+str(T10_14[i])+','+str(T15_19[i])+','\
            +str(T20_24[i])+','+str(T25_29[i])+','+str(T30_34[i])+','+str(T35_39[i])+','\
            +str(T40_44[i])+','+str(T45_49[i])+','+str(T50_54[i])+','+str(T55_59[i])+','\
            +str(T60_64[i])+','+str(T65_69[i])+','+str(T70_74[i])+','+str(T80_84[i])+','\
            +str(T90_94[i])+','+str(T95_99[i])
    outfile.write(Lifetable)
    outfile.write('\n')
for i in range(project_years):
    outfile2.write('year,age_group,male,female')
    outfile2.write('\n')
    for j in range(pk):
        Lifetable=str(i)+','+str(j)+','+str(male[i][j])+','+str(female[i][j])
        outfile2.write(Lifetable)
        outfile2.write('\n')
    
outfile.close()
outfile2.close()
infile.close()




      
