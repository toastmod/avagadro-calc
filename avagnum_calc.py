#!/usr/bin/env python3.7
from periodictable import elements

#generate table
print("generating table...")
ptable = {}
for e in elements:
    ptable[e.symbol] = {"name": e.name, "mass": e.mass, "num": e.number, "ions": e.ion}
print("done.\n")

# CUT HERE

avagnum = 6.022*(10**23)

def detectpwr(rcv):
    snd = 0.0 #final result
    cmd = ""
    pwr = 0.0
    
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

    if choice == 4:
        #moles to atoms
        moles = detectpwr(input("moles: "))
        print("atoms: "+str(moles_to_atoms(moles)))

    if choice == 5:
        #grams to moles
        grams = detectpwr(input("grams: "))
        molmass = detectpwr(input("molar mass: "))
        print("moles: "+str(grams/molmass))

    if choice == 6:
        #moles to grams
        moles = detectpwr(input("moles: "))
        molmass = detectpwr(input("molar mass: "))
        print("grams: "+str(moles*molmass))

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

    if choice == 8:
        fillorder = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p"]
        fillamt = {"s":2,"p":6,"d":10,"f":14}
        res = []
        atomnum = detectpwr(input("atomic number: "))
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
        print("electron config:")
        for i in res:
            print(i,end=" ")
        print("\n")
        print("MASTERINGCHEM friendly text:")
        for i in res:
            print(i,end="")
        print("\n")




    #========================================================

    if choice == 100:
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
        print("[enter]")

    #=================================================

    if choice == 0:
        #command test
        print(detectpwr(input("[$]>>")))
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




