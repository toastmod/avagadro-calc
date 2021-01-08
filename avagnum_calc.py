from periodictable import elements
import math
from decimal import Decimal, getcontext

tmp2 = 0

#generate table
print("generating table...")
ptable = {"STP":{"Pkpa":101.321,"Patm":1,"Tk":273.15,"Tf":32,"Tc":0,"molV":22.4,"molV1bar":22.7}}
for e in elements:
    ptable[e.symbol] = {"name": e.name, "mass": e.mass, "num": e.number, "ions": e.ion}
print("done.\n")

# CUT HERE

avagnum = 6.022*(10**23)

def log(x):
    
    return math.log(x,10)

def detectpwr(rcv):
    snd = 0.0 #final result
    cmd = ""
    pwr = 0.0

    if (str(rcv) == 'None') or (str(rcv) == ""):
        return ""
    
    #if command or not
    if str(rcv).find(":") >= 0:
        cmd = rcv.split(":")
        tmp = ptable[cmd[0]][cmd[1]]
        #if got number
        if str(tmp).replace(".","").isnumeric():
            if cmd[-1] == "e":
                print("yoo")
                pwr = float(input("*10^"))
                snd = float(tmp)*(10**pwr)
            else:
                snd = float(tmp)
        else:
            snd = tmp
    else:
        #no colon
        tmp = ""
        if str(rcv)[-1] == "e":
            tmp = str(rcv).replace("e","")
            #check if number (excluding decimal point)
            if str(tmp).replace(".","").isnumeric():
                pwr = float(input("*10^"))
                snd = float(tmp)*(10**pwr)
            else:
                #passthrough
                snd = rcv

        else:
            if str(rcv).replace(".","").isnumeric():
                #ensure number is output as float
                snd = float(rcv)
            else:
                #ensure word is output as str
                snd = str(rcv)


    return snd

def grams_to_atoms(grams,molmass):
    return (grams/molmass)*avagnum

def atoms_to_grams(atoms,molmass):
    return (atoms/avagnum)*molmass

def atoms_to_moles(atoms):
    return atoms/avagnum

def moles_to_atoms(moles):
    return moles*avagnum

def gen_e_config(atomnum):
        fillorder = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p"]
        fillamt = {"s":2,"p":6,"d":10,"f":14}
        res = []
        for fill in fillorder:
            if atomnum != 0:
                amt = fillamt[fill[1]]
                if amt > atomnum:
                    tmp = 0
                    while atomnum != 0:
                        tmp += 1
                        atomnum -= 1
                    res.append(fill+"^"+str(tmp))
                else:
                    atomnum -= amt
                    res.append(fill+"^"+str(amt))
        return res

def gen_qn(electron):
    global tmp2
    shell_vals = {"s":0,"p":1,"d":2,"f":3}
    level = int(electron[0])
    shell = electron[1]
    e_amt = int(electron[3])

    n = level
    l = shell_vals[shell]

    ML_fill = []
    tmp = -l
    if l != 0:
        while (tmp != l):
            ML_fill.append(tmp)
            tmp += 1
    else:
        ML_fill.append(0)

        tmp2 = 0
        down_spin = False
        for i in range(0,e_amt):
            if tmp2 >= len(ML_fill):
                tmp2 = 0
                down_spin = not down_spin
            else:
                tmp2 += 1

    ml = ML_fill[tmp2]
    ms = str("-"*(down_spin))+str("1/2")
    return [int(n),int(l),int(ml),ms]

# PV=nRT functions
# ideal gas constant
R = 0.0821

def pvnrt_v(P,T,n):
    return (n*R*T)/P

def pvnrt_p(V,n,T):
    return (n*R*T)/V

def pvnrt_n(P,V,T):
    return (P*V)/(R*T)

def pvnrt_t(P,V,n):
    (P*V)/(n*R)

#====================================================================

while 1==1:
    print("==[choose calc]====================")
    print("0. command output")
    print("1. grams -> atoms")
    print("2. atoms -> grams")
    print("3. atoms -> moles")
    print("4. moles -> atoms")
    print("5. grams -> moles")
    print("6. moles -> grams")
    print("7. quantum number validator")
    print("8. get electron configuration")
    print("9. get quantum number")
    print("10. get electrons and orbitals in shell")
    print("11. PVNRT")
    print("12. (LR_ratio,R_Liters) -> (L_moles,L_Grams)")
    print("13. (LR_ratio,R_moles)  -> (R_moles,R_Liters)")
    print("14. C1*V1=C2*V2")
    print("15. Beer's Law [ A=mc ]")
    print("16. Beer's Law [ A=Ebc ]")
    print("17. Beer's Law [ A=log(T) | A1/A2=C1/C2 ]")
    print("18. Determine the Mole Ratio (from moles)")
    print("19. Determine the Mole Ratio (from grams/molmass)")
    print("20. mmH2O -> mmHg")
    print("21. mmHg -> mmH2O")
    print("22. rate of reaction (given a timetable)")
    print()
    choice = int(input("[#]: "))
    print("===================================")
    if choice == 1:
        #grams to atoms
        grams = detectpwr(input("grams: "))
        molmass = detectpwr(input("molar mass: "))
        atoms = grams_to_atoms(grams,molmass)
        print("atoms: "+str(atoms))
        input("[enter]")

    if choice == 2:
        #atoms to grams
        atoms = detectpwr(input("atoms: "))
        molmass = detectpwr(input("molar mass: "))
        grams = atoms_to_grams(atoms,molmass)
        print("grams: "+str(grams))
        input("[enter]")

    if choice == 3:
        #atoms to moles
        atoms = detectpwr(input("atoms: "))
        print("moles: "+str(atoms_to_moles(atoms)))
        input("[enter]")

    if choice == 4:
        #moles to atoms
        moles = detectpwr(input("moles: "))
        print("atoms: "+str(moles_to_atoms(moles)))
        input("[enter]")

    if choice == 5:
        #grams to moles
        grams = detectpwr(input("grams: "))
        molmass = detectpwr(input("molar mass: "))
        print("moles: "+str(grams/molmass))
        input("[enter]")

    if choice == 6:
        #moles to grams
        moles = detectpwr(input("moles: "))
        molmass = detectpwr(input("molar mass: "))
        print("grams: "+str(moles*molmass))
        input("[enter]")

    if choice == 7:
        #quantum number validator
        entr = input("enter qn (list with commas): ")
        qn = entr.split(",")

        n = int(qn[0])
        l = int(qn[1])
        ml = int(qn[2])
        ms = qn[3]

        flag = False
        #print(n >= 1)
        if n >= 1:
            #print((l > 0),(l < (n-1)))
            if (l >= 0) and (l <= (n-1)):
                #print((ml <= l),(ml >= -l))
                if (ml <= l) and (ml >= -l):
                    #print((ms == "1/2"),(ms == "-1/2"))
                    if (ms == "1/2") or (ms == "-1/2"):
                        flag = True
        if flag==False:
            print("QN INVALID")
        else:   
            print("QN VALID :)")
        input("[enter]")

    if choice == 8:
        res = gen_e_config(detectpwr(input("atomic number: ")))
        fillamt = {"s":2,"p":6,"d":10,"f":14}

        magnetism = "diamagnetic"
        tmp = res[-1].split("^")
        if int(tmp[1]) < fillamt[tmp[0][1]]:
            magnetism = "paramagnetic"


        print("electron config:")
        for i in res:
            print(i,end=" ")
        print("\n")
        print("MASTERINGCHEM friendly text:")
        for i in res:
            print(i,end="")
        print("\n\nelement is "+magnetism+"\n")
        input("[enter]")

    if choice == 9:
        electrons = gen_e_config(detectpwr(input("atomic number: ")))
        print("POSSIBLE QUANTUM NUMBERS:")
        print("")
        for electron in electrons:
            qn = gen_qn(electron)
            print(electron+": "+str(qn).replace("[","").replace("]","").replace("'",""))
            if qn[-1][0] == "-":
                qn2 = qn
                qn2[-1] = qn2[-1].replace("-","")
                print("      "+str(qn2).replace("[","").replace("]","").replace("'",""))
            print("")
        input("[enter]")

    if choice == 10:
        subshells = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p"]
        eamt = {"s":2,"p":6,"d":10,"f":14}
        tmp = 0
        inp = input("n = ")
        for s in subshells:
            if inp == s[0]:
                tmp += eamt[s[1]]
        print(str(int(tmp))+" electrons")
        print(str(int(tmp/2))+" orbitals")
        input("[enter]")

    if choice == 11:
        P = detectpwr(input("P="))
        V = detectpwr(input("V="))
        n = detectpwr(input("n="))
        T = detectpwr(input("T="))

        if P == "":
            P = pvnrt_p(float(V),float(n),float(T))
            print("P -> "+str(P))
        if V == "":
            V = pvnrt_v(float(P),float(n),float(T))
            print("V -> "+str(V))
        if n == "":
            n = pvnrt_n(float(P),float(V),float(T))
            print("n -> "+str(n))
        if T == "":
            T = pvnrt_t(float(P),float(V),float(n))
            print("T -> "+str(T))

    
    if choice==12:
        Lr = detectpwr(input("Left ratio: "))
        Rr = detectpwr(input("Right ratio: "))
        V = detectpwr(input("Right Volume: "))
        Rm = V*1000
        print("PVNRT environment:\n")
        P = detectpwr(input("P="))
        n = ""
        T = detectpwr(input("T="))

        if P == "":
            P = pvnrt_p(float(V),float(n),float(T))
            print("P -> "+str(P))
        if V == "":
            V = pvnrt_v(float(P),float(n),float(T))
            print("V -> "+str(V))
        if n == "":
            n = pvnrt_n(float(P),float(V),float(T))
            print("n -> "+str(n))
        if T == "":
            T = pvnrt_t(float(P),float(V),float(n))
            print("T -> "+str(T))

        answer_mol = (n*Lr)/Rr
        answer_grams = atoms_to_grams(moles_to_atoms(answer_mol),detectpwr(input("Left molar mass: ")))
        print("Left:")
        print("\tMoles> "+str(answer_mol))
        print("\tGrams> "+str(answer_grams))

    if choice==13:
        Lr = detectpwr(input("Left ratio: "))
        Rr = detectpwr(input("Right ratio: "))
        Ln = detectpwr(input("Left moles used: "))
        # get right ratio
        n = (Ln*Rr)/Lr

        print("PVNRT environment:\n")
        P = detectpwr(input("P="))
        V = ""
        T = detectpwr(input("T="))

        if P == "":
            P = pvnrt_p(float(V),float(n),float(T))
            print("P -> "+str(P))
        if V == "":
            V = pvnrt_v(float(P),float(n),float(T))
            print("V -> "+str(V))
        if n == "":
            n = pvnrt_n(float(P),float(V),float(T))
            print("n -> "+str(n))
        if T == "":
            T = pvnrt_t(float(P),float(V),float(n))
            print("T -> "+str(T))

        answer_mol = n
        #answer_liters = atoms_to_grams(moles_to_atoms(answer_mol),detectpwr(input("Right molar mass: ")))/1000
        print("Right:")
        print("\tMoles> "+str(answer_mol))
        print("\tLiters> "+str(V))

    if choice==14:
        #c1*v1=c2*v2
        c1 = detectpwr(input("C1="))
        v1 = detectpwr(input("V1="))
        c2 = detectpwr(input("C2="))
        v2 = detectpwr(input("V2="))

        if c1=="":
            print("\nC1 -> "+str((c2*v2)/v1))

        if c2=="":
            print("\nC2 -> "+str((c1*v1)/v2))

        if v1=="":
            print("\nC1 -> "+str((c2*v2)/c1))

        if v2=="":
            print("\nC2 -> "+str((c1*v1)/c2))

    if choice==15:
        A = detectpwr(input("A=")) # absorvity
        c = detectpwr(input("c[concentration]="))
        m = detectpwr(input("m (slope)="))

        if A=="":
            A=m*c
            print("A -> "+str(A))

        if c=="":
            c=A/m
            print("c -> "+str(c))

    if choice==16:
        A = detectpwr(input("A=")) # absorvity
        c = detectpwr(input("c[concentration]="))
        E = detectpwr(input("E="))
        b = detectpwr(input("b[molar absorp.]="))

        if A=="":
            A=E*b*c
            print("A -> "+str(A))

        if E=="":
            E=A/(b*c)
            print("E -> "+str(E))

        if b=="":
            b=A/(E*c)
            print("b -> "+str(b))

        if c=="":
            c=A/(E*b)
            print("c -> "+str(c))

    if choice==17:
        
        a1 = detectpwr(input("Absorbance    of [1]="))
        c1 = detectpwr(input("Concentration of [1]="))
        t1 = detectpwr(input("Tranmittance% of [1]="))

        a2 = detectpwr(input("Absorbance    of [2]="))
        c2 = detectpwr(input("Concentration of [2]="))
        t2 = detectpwr(input("Tranmittance% of [2]="))


        if a1=="":
            if (t1!=""):
                a1 = log(100/t1)
                print("A1 -> "+str(a1))

            elif (a2 == ""):
                # caculate a2 with t2
                a2 = log(100/t2)
                print("A2 -> "+str(a2))
                a1 = (c1*a2)/c2
                print("A1 -> "+str(a1))

        if a2=="":
            if (t2!=""):
                a2 = log(100/t2)
                print("A2 -> "+str(a2))

            elif (a1 == ""):
                # caculate a1 with t1
                a1 = log(100/t1)
                print("A1 -> "+str(a1))
                a2 = (c2*a1)/c1
                print("A2 -> "+str(a2))

        if c1=="":
            if (a2 == ""):
                # caculate a2 with t2
                a2 = log(100/t2)
                print("A2 -> "+str(a2))

            elif (a1 == ""):
                # caculate a1 with t1
                a1 = log(100/t1)
                print("A1 -> "+str(a1))

            c1 = (a1*c2)/a2
            print("C1 -> "+str(c1))

        if c2=="":
            if (a2 == ""):
                # caculate a2 with t2
                a2 = log(100/t2)
                print("A2 -> "+str(a2))

            elif (a1 == ""):
                # caculate a1 with t1
                a1 = log(100/t1)
                print("A1 -> "+str(a1))

            c2 = (a2*c1)/a1
            print("C2 -> "+str(c2))

    if choice==18:
        n1 = detectpwr(input("left [moles]:"))
        n2 = detectpwr(input("right [moles]:"))

        lowest = 0
        highest = 0

        if n1 > n2:
            lowest = n2
            highest = n1
            lowscale = 1/lowest
            r1 = highest*lowscale
            r2 = lowest*lowscale
        else:
            lowest = n1
            highest = n2
            lowscale = 1/lowest
            r2 = highest*lowscale
            r1 = lowest*lowscale

        print("Ratio (exact): "+str(r1)+":"+str(r2))
        print("Ratio (trunc): "+str(int(r1))+":"+str(int(r2)))

    if choice==19:
        m1 = detectpwr(input("left [grams]:"))
        molmass1 = detectpwr(input("molar mass: "))
        m2 = detectpwr(input("right [grams]:"))
        molmass2 = detectpwr(input("molar mass: "))

        n1 = atoms_to_moles((grams_to_atoms(m1,molmass1)))
        n2 = atoms_to_moles((grams_to_atoms(m2,molmass2)))

        lowest = 0
        highest = 0

        if n1 > n2:
            lowest = n2
            highest = n1
            lowscale = 1/lowest
            r1 = highest*lowscale
            r2 = lowest*lowscale
        else:
            lowest = n1
            highest = n2
            lowscale = 1/lowest
            r2 = highest*lowscale
            r1 = lowest*lowscale

        print("Ratio (exact): "+str(r1)+":"+str(r2))
        print("Ratio (trunc): "+str(int(r1))+":"+str(int(r2)))

    if choice==20:
        bar = detectpwr(input("bar (optional): "))
        mmh2o = detectpwr(input("mmH2O: "))
        print("mmHg -> "+str(bar-(0.0735559*mmh2o)))

    if choice==21:
        bar = detectpwr(input("bar (optional): "))
        mmhg = bar-detectpwr(input("mmHg: "))
        print("mmH2O -> "+str((mmhg/0.0735559)))

    if choice==22:
        x = detectpwr(input("Formation@timeX:"))
        y = detectpwr(input("Formation@timey:"))
        print("change = "+str(abs(x-y)))
        
        
        




    #========================================================

    if choice == 9999:
        #simplify mole ratio
        elem_mol = {} #collect element and moles here
        highest_tenth = 1;
        tmp = 0;

        #kickoff
        name = detectpwr(input("enter element name [or blank to stop]: "))
        while name != "":
            mols = detectpwr(input("enter grams per mole of element      : "))

            #check how high the decimal needs to be trunc'd
            tmp_mols = mols
            nodec = float(str(tmp_mols).replace(".",""))
            index = 0; #start at 1 to avoid div by 0
            loop_flag = True
            #there are way more efficient ways to do this but whatever
            while loop_flag:
                print(tmp_mols)
                if tmp_mols != nodec:
                    tmp_mols *= 10
                    index += 1
                if tmp_mols == nodec:
                    if index > highest_tenth:
                        highest_tenth = index
                    loop_flag = False

            #write to dict
            elem_gpm[name] = mols
            #next element (+ test for input continuation)
            name = detectpwr(input("enter element name [or blank to stop]: "))

        #once collection is done transform to a list for indexing
        elem_mol = list(elem_mol.items())

        # multiply all to the highest tenth
        for i in range(0,len(elem_mol)):
            # turn tuple into a list first
            elem_mol[i] = list(elem_mol[i])
            elem_mol[i][1] *= (10**highest_tenth)

        # now brute force the common factor
        loop_flag2 = True
        mod_flag = True #true while there is modulo
        index2 = 1;
        while loop_flag2:
            # check if mod flag holds up
            mod_flag = False
            for e in elem_mol:
                if (e[1]%index2) != 0:
                    mod_flag = True 

            #if modflag is not held up
            if mod_flag == False:
                #we got em
                loop_flag2 = False
        
        # on loop break
        print("common factor: "+str(index2))
        for e in elem_mol:
            print(str(e[0])+str(e[1]/index2)+":",end="")

        parts = len(elem_mol)
        input("[enter]")

    #=================================================

    if choice == 0:
        #command test
        inp = input("[$]>>")
        inp2 = "0"
        sum_this = []
        sum = 0.0
        if inp == "add":
            while inp2 != "":
                sum_this.append(inp2)
                inp2 = str(detectpwr(e))
            for e in sum_this:
                sum += float(e)
            print(sum)
        else:
            print(detectpwr(inp))
        input("[enter]")

    if choice == 999:
        #UNFINISHED
        elem_gpm = {}
        name = "PLACEHOLDER"
        gpm = 0
        while name != "":
            name = detectpwr(input("enter element name [or blank to stop]: "))
            gpm = detectpwr(input("enter grams per mole of element      : "))
            elem_gpm[name] = gpm
        parts = len(elem_gpm)
        print("")




