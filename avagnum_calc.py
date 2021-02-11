import math
import json
from decimal import Decimal, getcontext

tmp2 = 0

time0 = 0

# generate table
from periodictable import elements

ptable = {"STP": {"Pkpa": 101.321, "Patm": 1, "Tk": 273.15, "Tf": 32, "Tc": 0, "molV": 22.4, "molV1bar": 22.7}}
for e in elements:
    ptable[e.symbol] = {"name": e.name, "mass": e.mass, "num": e.number, "ions": tuple(e.ions), "isotopes": list(e.isotopes), "charge": e.charge}
# print("loading data...")
#dataf = open(".\\data.txt","r")
#datar = dataf.readline()
#print(datar)
#ptable = json.loads(datar)


avagnum = 6.022 * (10 ** 23)


def log(x):
    return math.log(x, 10)


def ln(x):
    return math.log(x)


def un_log(y):
    return 10 ** y


def un_ln(y):
    return math.exp(y)


def detectpwr(rcv):
    snd = 0.0  # final result
    cmd = ""
    pwr = 0.0

    if str(rcv) == "eq":
        return eval(input("[eq]>> "))


    if (str(rcv) == 'None') or (str(rcv) == ""):
        return ""

    # if command or not
    if str(rcv).find(":") >= 0:
        cmd = rcv.split(":")
        tmp = ptable[cmd[0]][cmd[1]]
        # if got number
        if str(tmp).replace(".", "").isnumeric():
            if cmd[-1] == "e":
                pwr = float(input("*10^"))
                snd = float(tmp) * (10 ** pwr)
            else:
                snd = float(tmp)
        else:
            snd = tmp
    else:
        # no colon
        snd = ""
        if str(rcv)[-1] == "e":
            tmp = str(rcv).replace("e", "")
            # check if number (excluding decimal point)
            if str(tmp).replace(".", "").isnumeric():
                pwr = float(input("*10^"))
                snd = float(tmp) * (10 ** pwr)
                return snd
            
            else:
                # passthrough
                snd = rcv
        elif str(rcv)[-1] == "C":
            print("[Celsius to Kelvin]")
            tmp = str(rcv).replace("C", "")
            return float(float(tmp) + 273.15)
            print(snd)
        elif str(rcv)[-1] == "K":
            print("[Kelvin to Celsius]")
            tmp = str(rcv).replace("K", "")
            return float(float(tmp) - 273.15)

        else:
            if str(rcv).replace(".", "").isnumeric():
                # ensure number is output as float
                snd = float(rcv)
            else:
                # ensure word is output as str
                snd = str(rcv)

    return snd


def grams_to_atoms(grams, molmass):
    return (grams / molmass) * avagnum


def atoms_to_grams(atoms, molmass):
    return (atoms / avagnum) * molmass


def atoms_to_moles(atoms):
    return atoms / avagnum


def moles_to_atoms(moles):
    return moles * avagnum


def gen_e_config(atomnum):
    fillorder = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f",
                 "6d", "7p"]
    fillamt = {"s": 2, "p": 6, "d": 10, "f": 14}
    res = []
    for fill in fillorder:
        if atomnum != 0:
            amt = fillamt[fill[1]]
            if amt > atomnum:
                tmp = 0
                while atomnum != 0:
                    tmp += 1
                    atomnum -= 1
                res.append(fill + "^" + str(tmp))
            else:
                atomnum -= amt
                res.append(fill + "^" + str(amt))
    return res


def gen_qn(electron):
    global tmp2
    shell_vals = {"s": 0, "p": 1, "d": 2, "f": 3}
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
        for i in range(0, e_amt):
            if tmp2 >= len(ML_fill):
                tmp2 = 0
                down_spin = not down_spin
            else:
                tmp2 += 1

    ml = ML_fill[tmp2]
    ms = str("-" * (down_spin)) + str("1/2")
    return [int(n), int(l), int(ml), ms]


# PV=nRT functions
# ideal gas constant
R = 0.0821


def pvnrt_v(P, T, n):
    return (n * R * T) / P


def pvnrt_p(V, n, T):
    return (n * R * T) / V


def pvnrt_n(P, V, T):
    return (P * V) / (R * T)


def pvnrt_t(P, V, n):
    (P * V) / (n * R)


# Nth-Order Reaction

def zero_order():
    A0 = detectpwr(input("[A]0="))
    k = detectpwr(input("k="))
    t = detectpwr(input("t="))

    tmp_half = input("[A]t ( type 'half:#' to use (1/#)*[A]0 ) = ").split(":")

    if tmp_half[0] == "half":
        At = A0 / float(tmp_half[1])
    else:
        At = detectpwr(tmp_half[0])

    At = detectpwr()

    if At == "":
        t = t - time0
        At = ((-k * t) + (A0))
        print("\n[A]t = " + str(At))
    if A0 == "":
        t = t - time0
        A0 = ((At) + (k * t))
        print("\n[A]0 = " + str(A0))
    if k == "":
        t = t - time0
        k = -((At) - (A0)) / t
        print("\nk = " + str(k))
    if t == "":
        t = ((At) - (A0)) / -k
        t = t - time0
        print("\nt = " + str(t))


def first_order():
    A0 = detectpwr(input("[A]0="))
    k = detectpwr(input("k="))
    t = detectpwr(input("t="))
    At = detectpwr(input("[A]t="))

    if At == "":
        t = t - time0
        A0 = ln(A0)
        At = ((-k * t) + (A0))
        At = un_ln(At)
        print("\n[A]t = " + str(At))
    if A0 == "":
        t = t - time0
        At = ln(At)
        A0 = ((At) + (k * t))
        A0 = un_ln(A0)
        print("\n[A]0 = " + str(A0))
    if k == "":
        t = t - time0
        At = ln(At)
        A0 = ln(A0)
        k = -((At) - (A0)) / t
        print("\nk = " + str(k))
    if t == "":
        At = ln(At)
        A0 = ln(A0)
        t = -((At) - (A0)) / k
        t = t - time0
        print("\nt = " + str(t))


def second_order():
    A0 = detectpwr(input("[A]0="))
    k = detectpwr(input("k="))
    t = detectpwr(input("t="))
    At = detectpwr(input("[A]t="))

    # At = 1/((k*t)+(1/A0))

    if At == "":
        t = t - time0
        At = 1 / ((k * t) + (1 / A0))
        print("\n[A]t = " + str(At))
    if A0 == "":
        t = t - time0
        A0 = 1 / ((1 / At) - (k * t))
        print("\n[A]0 = " + str(A0))
    if k == "":
        t = t - time0
        k = ((1 / At) - (1 / A0)) / t
        print("\nk = " + str(k))
    if t == "":
        t = ((1 / At) - (1 / A0)) / k
        t = t - time0
        print("\nt = " + str(t))


# ====================================================================

while 1 == 1:
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
    print("22. rate of formation/decomp (given time+formation)")
    print("23. Formation Rate -> Consumption Rate (given moles)")
    print("24. Nth-order reaction formula")
    print("25. Half-life formula")
    print("26. Arrhenius Equation")
    print("27. Kc <-> Kp (unfinished)")
    print("28. Kw = [H+]*[OH-]")
    
    # print("24. Rate Law")
    print()
    choice = int(input("[#]: "))
    print("===================================")
    if choice == 1:
        # grams to atoms
        grams = detectpwr(input("grams: "))
        molmass = detectpwr(input("molar mass: "))
        atoms = grams_to_atoms(grams, molmass)
        print("atoms: " + str(atoms))
        input("[enter]")

    if choice == 2:
        # atoms to grams
        atoms = detectpwr(input("atoms: "))
        molmass = detectpwr(input("molar mass: "))
        grams = atoms_to_grams(atoms, molmass)
        print("grams: " + str(grams))
        input("[enter]")

    if choice == 3:
        # atoms to moles
        atoms = detectpwr(input("atoms: "))
        print("moles: " + str(atoms_to_moles(atoms)))
        input("[enter]")

    if choice == 4:
        # moles to atoms
        moles = detectpwr(input("moles: "))
        print("atoms: " + str(moles_to_atoms(moles)))
        input("[enter]")

    if choice == 5:
        # grams to moles
        grams = detectpwr(input("grams: "))
        molmass = detectpwr(input("molar mass: "))
        print("moles: " + str(grams / molmass))
        input("[enter]")

    if choice == 6:
        # moles to grams
        moles = detectpwr(input("moles: "))
        molmass = detectpwr(input("molar mass: "))
        print("grams: " + str(moles * molmass))
        input("[enter]")

    if choice == 7:
        # quantum number validator
        entr = input("enter qn (list with commas): ")
        qn = entr.split(",")

        n = int(qn[0])
        l = int(qn[1])
        ml = int(qn[2])
        ms = qn[3]

        flag = False
        # print(n >= 1)
        if n >= 1:
            # print((l > 0),(l < (n-1)))
            if (l >= 0) and (l <= (n - 1)):
                # print((ml <= l),(ml >= -l))
                if (ml <= l) and (ml >= -l):
                    # print((ms == "1/2"),(ms == "-1/2"))
                    if (ms == "1/2") or (ms == "-1/2"):
                        flag = True
        if flag == False:
            print("QN INVALID")
        else:
            print("QN VALID :)")
        input("[enter]")

    if choice == 8:
        res = gen_e_config(detectpwr(input("atomic number: ")))
        fillamt = {"s": 2, "p": 6, "d": 10, "f": 14}

        magnetism = "diamagnetic"
        tmp = res[-1].split("^")
        if int(tmp[1]) < fillamt[tmp[0][1]]:
            magnetism = "paramagnetic"

        print("electron config:")
        for i in res:
            print(i, end=" ")
        print("\n")
        print("MASTERINGCHEM friendly text:")
        for i in res:
            print(i, end="")
        print("\n\nelement is " + magnetism + "\n")
        input("[enter]")

    if choice == 9:
        electrons = gen_e_config(detectpwr(input("atomic number: ")))
        print("POSSIBLE QUANTUM NUMBERS:")
        print("")
        for electron in electrons:
            qn = gen_qn(electron)
            print(electron + ": " + str(qn).replace("[", "").replace("]", "").replace("'", ""))
            if qn[-1][0] == "-":
                qn2 = qn
                qn2[-1] = qn2[-1].replace("-", "")
                print("      " + str(qn2).replace("[", "").replace("]", "").replace("'", ""))
            print("")
        input("[enter]")

    if choice == 10:
        subshells = ["1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s",
                     "5f", "6d", "7p"]
        eamt = {"s": 2, "p": 6, "d": 10, "f": 14}
        tmp = 0
        inp = input("n = ")
        for s in subshells:
            if inp == s[0]:
                tmp += eamt[s[1]]
        print(str(int(tmp)) + " electrons")
        print(str(int(tmp / 2)) + " orbitals")
        input("[enter]")

    if choice == 11:
        P = detectpwr(input("P="))
        V = detectpwr(input("V="))
        n = detectpwr(input("n="))
        T = detectpwr(input("T="))

        if P == "":
            P = pvnrt_p(float(V), float(n), float(T))
            print("P -> " + str(P))
        if V == "":
            V = pvnrt_v(float(P), float(n), float(T))
            print("V -> " + str(V))
        if n == "":
            n = pvnrt_n(float(P), float(V), float(T))
            print("n -> " + str(n))
        if T == "":
            T = pvnrt_t(float(P), float(V), float(n))
            print("T -> " + str(T))
        input("[enter]")

    if choice == 12:
        Lr = detectpwr(input("Left ratio: "))
        Rr = detectpwr(input("Right ratio: "))
        V = detectpwr(input("Right Volume: "))
        Rm = V * 1000
        print("PVNRT environment:\n")
        P = detectpwr(input("P="))
        n = ""
        T = detectpwr(input("T="))

        if P == "":
            P = pvnrt_p(float(V), float(n), float(T))
            print("P -> " + str(P))
        if V == "":
            V = pvnrt_v(float(P), float(n), float(T))
            print("V -> " + str(V))
        if n == "":
            n = pvnrt_n(float(P), float(V), float(T))
            print("n -> " + str(n))
        if T == "":
            T = pvnrt_t(float(P), float(V), float(n))
            print("T -> " + str(T))

        answer_mol = (n * Lr) / Rr
        answer_grams = atoms_to_grams(moles_to_atoms(answer_mol), detectpwr(input("Left molar mass: ")))
        print("Left:")
        print("\tMoles> " + str(answer_mol))
        print("\tGrams> " + str(answer_grams))
        input("[enter]")

    if choice == 13:
        Lr = detectpwr(input("Left ratio: "))
        Rr = detectpwr(input("Right ratio: "))
        Ln = detectpwr(input("Left moles used: "))
        # get right ratio
        n = (Ln * Rr) / Lr

        print("PVNRT environment:\n")
        P = detectpwr(input("P="))
        V = ""
        T = detectpwr(input("T="))

        if P == "":
            P = pvnrt_p(float(V), float(n), float(T))
            print("P -> " + str(P))
        if V == "":
            V = pvnrt_v(float(P), float(n), float(T))
            print("V -> " + str(V))
        if n == "":
            n = pvnrt_n(float(P), float(V), float(T))
            print("n -> " + str(n))
        if T == "":
            T = pvnrt_t(float(P), float(V), float(n))
            print("T -> " + str(T))

        answer_mol = n
        # answer_liters = atoms_to_grams(moles_to_atoms(answer_mol),detectpwr(input("Right molar mass: ")))/1000
        print("Right:")
        print("\tMoles> " + str(answer_mol))
        print("\tLiters> " + str(V))
        input("[enter]")

    if choice == 14:
        # c1*v1=c2*v2
        c1 = detectpwr(input("C1="))
        v1 = detectpwr(input("V1="))
        c2 = detectpwr(input("C2="))
        v2 = detectpwr(input("V2="))

        if c1 == "":
            print("\nC1 -> " + str((c2 * v2) / v1))

        if c2 == "":
            print("\nC2 -> " + str((c1 * v1) / v2))

        if v1 == "":
            print("\nC1 -> " + str((c2 * v2) / c1))

        if v2 == "":
            print("\nC2 -> " + str((c1 * v1) / c2))
        input("[enter]")

    if choice == 15:
        A = detectpwr(input("A="))  # absorvity
        c = detectpwr(input("c[concentration]="))
        m = detectpwr(input("m (slope)="))

        if A == "":
            A = m * c
            print("A -> " + str(A))

        if c == "":
            c = A / m
            print("c -> " + str(c))
        input("[enter]")

    if choice == 16:
        A = detectpwr(input("A="))  # absorvity
        c = detectpwr(input("c[concentration]="))
        E = detectpwr(input("E="))
        b = detectpwr(input("b[molar absorp.]="))

        if A == "":
            A = E * b * c
            print("A -> " + str(A))

        if E == "":
            E = A / (b * c)
            print("E -> " + str(E))

        if b == "":
            b = A / (E * c)
            print("b -> " + str(b))

        if c == "":
            c = A / (E * b)
            print("c -> " + str(c))
        input("[enter]")

    if choice == 17:

        a1 = detectpwr(input("Absorbance    of [1]="))
        c1 = detectpwr(input("Concentration of [1]="))
        t1 = detectpwr(input("Tranmittance% of [1]="))

        a2 = detectpwr(input("Absorbance    of [2]="))
        c2 = detectpwr(input("Concentration of [2]="))
        t2 = detectpwr(input("Tranmittance% of [2]="))

        if a1 == "":
            if (t1 != ""):
                a1 = log(100 / t1)
                print("A1 -> " + str(a1))

            elif (a2 == ""):
                # caculate a2 with t2
                a2 = log(100 / t2)
                print("A2 -> " + str(a2))
                a1 = (c1 * a2) / c2
                print("A1 -> " + str(a1))

        if a2 == "":
            if (t2 != ""):
                a2 = log(100 / t2)
                print("A2 -> " + str(a2))

            elif (a1 == ""):
                # caculate a1 with t1
                a1 = log(100 / t1)
                print("A1 -> " + str(a1))
                a2 = (c2 * a1) / c1
                print("A2 -> " + str(a2))

        if c1 == "":
            if (a2 == ""):
                # caculate a2 with t2
                a2 = log(100 / t2)
                print("A2 -> " + str(a2))

            elif (a1 == ""):
                # caculate a1 with t1
                a1 = log(100 / t1)
                print("A1 -> " + str(a1))

            c1 = (a1 * c2) / a2
            print("C1 -> " + str(c1))

        if c2 == "":
            if (a2 == ""):
                # caculate a2 with t2
                a2 = log(100 / t2)
                print("A2 -> " + str(a2))

            elif (a1 == ""):
                # caculate a1 with t1
                a1 = log(100 / t1)
                print("A1 -> " + str(a1))

            c2 = (a2 * c1) / a1
            print("C2 -> " + str(c2))
        input("[enter]")

    if choice == 18:
        n1 = detectpwr(input("left [moles]:"))
        n2 = detectpwr(input("right [moles]:"))

        lowest = 0
        highest = 0

        if n1 > n2:
            lowest = n2
            highest = n1
            lowscale = 1 / lowest
            r1 = highest * lowscale
            r2 = lowest * lowscale
        else:
            lowest = n1
            highest = n2
            lowscale = 1 / lowest
            r2 = highest * lowscale
            r1 = lowest * lowscale

        print("Ratio (exact): " + str(r1) + ":" + str(r2))
        print("Ratio (trunc): " + str(int(r1)) + ":" + str(int(r2)))
        input("[enter]")

    if choice == 19:
        m1 = detectpwr(input("left [grams]:"))
        molmass1 = detectpwr(input("molar mass: "))
        m2 = detectpwr(input("right [grams]:"))
        molmass2 = detectpwr(input("molar mass: "))

        n1 = atoms_to_moles((grams_to_atoms(m1, molmass1)))
        n2 = atoms_to_moles((grams_to_atoms(m2, molmass2)))

        lowest = 0
        highest = 0

        if n1 > n2:
            lowest = n2
            highest = n1
            lowscale = 1 / lowest
            r1 = highest * lowscale
            r2 = lowest * lowscale
        else:
            lowest = n1
            highest = n2
            lowscale = 1 / lowest
            r2 = highest * lowscale
            r1 = lowest * lowscale

        print("Ratio (exact): " + str(r1) + ":" + str(r2))
        print("Ratio (trunc): " + str(int(r1)) + ":" + str(int(r2)))
        input("[enter]")

    if choice == 20:
        bar = detectpwr(input("bar (optional): "))
        mmh2o = detectpwr(input("mmH2O: "))
        print("mmHg -> " + str(bar - (0.0735559 * mmh2o)))
        input("[enter]")

    if choice == 21:
        bar = detectpwr(input("bar (optional): "))
        mmhg = bar - detectpwr(input("mmHg: "))
        print("mmH2O -> " + str((mmhg / 0.0735559)))
        input("[enter]")

    if choice == 22:
        floop = True
        ft = []
        while floop:
            form = detectpwr(input("formation: "))
            if form == "":
                floop = False
            time = detectpwr(input("time: "))
            if time == "":
                floop = False

            if floop == True:
                ft.append([form, time])

        print("*********")
        print("")
        sum_ft = 0
        for i in range(0, len(ft) - 1):
            delta_f = ft[i + 1][0] - ft[i][0]
            delta_t = ft[i + 1][1] - ft[i][1]
            delta_ft = delta_f / delta_t
            sum_ft += delta_ft
            print(str(delta_f / delta_t), end=",")

        print("")
        average = sum_ft / (len(ft) - 1)
        print("average=" + str(average))
        print("decomp.=" + str(2 * average))
        input("[enter]")

    if choice == 23:
        print("[EQ LEFT SIDE]")
        c_mol = detectpwr(input("Consumed moles: "))
        c_rate = detectpwr(input("Consumption rate: "))
        print("[EQ RIGHT SIDE]")
        f_mol = detectpwr(input("Formed moles: "))
        f_rate = detectpwr(input("Formation rate: "))
        print("************")
        if f_rate == "":
            f_rate = (c_rate * f_mol) / c_mol
            print("> Formation rate: " + str(f_rate))
        if c_rate == "":
            c_rate = (f_rate * c_mol) / f_mol
            print("> Consumption rate: " + str(c_rate))

        input("[enter]")

    if choice == 24:
        print("[0] zero order")
        print("[1] first order")
        print("[2] second order")
        choice2 = int(input("[#]: "))

        tmp0 = detectpwr(input("set initial time @ [A]0 (default=0):"))

        if tmp0 == "":
            time0 = 0
        else:
            time0 = tmp0

        if choice2 == 0:
            zero_order()
        if choice2 == 1:
            first_order()
        if choice2 == 2:
            second_order()

        input("[enter]")

    if choice == 25:
        # halflife kinetics
        c = 0.693
        thalf = detectpwr(input("t1/2="))
        k = detectpwr(input("k="))

        if thalf == "":
            thalf = c / k
            print("\nt1/2 = " + str(thalf))

        if k == "":
            k = 1 / (thalf / c)
            print("\nk = " + str(k))

        input("[enter]")

    if choice == 26:
        R = 8.314
        print("[1] Non-Linear")
        print("[2] ln(k2/k1)=(Ea/R)((1/T1)-(1/T2))")

        arrchoice = input("[#]>")

        if arrchoice == "1":

            Ea = detectpwr(input("Ea="))
            T = detectpwr(input("T="))
            A = detectpwr(input("A="))
            k = detectpwr(input("k="))

            print("************************")
            if Ea == "":
                Ea = ln(k / A) * R * T
                print("Ea=" + str(Ea))
            if T == "":
                T = Ea / (ln(k / A) * R)
                print("T=" + str(T))
            if k == "":
                k = A * un_ln(-Ea / (R * T))
                print("k=" + str(k))
            if A == "":
                A = k / un_ln(-Ea / (R * T))
                print("A=" + str(A))

        if (arrchoice == "2") or (arrchoice == "3"):

            # =========T and K entry==========
            print("==[T]==")
            print("[1] enter full ((1/T1)-(1/T2)) value")
            print("[2] enter T1 and T2 individually")
            tchoice = input("[#]>")
            if tchoice == "1":
                tinvdelta = detectpwr(input("((1/T1)-(1/T2))="))

            if tchoice == "2":
                T1 = detectpwr(input("T1="))
                T2 = detectpwr(input("T2="))

                if (T1 == "") or (T2 == ""):
                    tinvdelta = ""
                else:
                    tinvdelta = ((1 / T1) - (1 / T2))

            print("==[K]==")
            print("[1] enter full ln(k2/k1) value")
            print("[2] enter k1 and k2 individually")
            kchoice = input("[#]>")
            if kchoice == "1":
                lnk2k1 = detectpwr(input("ln(k2/k1)="))
            if kchoice == "2":
                k1 = detectpwr(input("k1="))
                k2 = detectpwr(input("k2="))
                if (k1 == "") or (k2 == ""):
                    lnk2k1 = ""
                else:
                    lnk2k1 = ln(k2/k1)
            if arrchoice == "2":
                Ea = detectpwr(input("Ea="))
                print(Ea)

            #        ==== EVAL ====
            print("**********************")
            # K
            if lnk2k1 == "":
                k2flag = True
                if k1 == "":
                    if k2 == "":
                        k2flag = False # k2 was already checked
                        # lnk2k1 is being looked for
                        lnk2k1 = (Ea/R)*tinvdelta
                        print("ln(k2/k1)=" + str(lnk2k1))
                    else:
                        # just k1 is being looked for
                        if arrchoice == "2":
                            k1 = 1/(un_ln((Ea/R)*tinvdelta)/k2)
                        if arrchoice == "3":
                            k1 = 1/(un_ln(tinvdelta)/k2)
                        print("k1=" + str(k1))


                if (k2 == "") and k2flag:
                    # just k2 is being looked for
                    if arrchoice == "2":
                        k2 = (un_ln((Ea / R) * tinvdelta) * k1)
                    if arrchoice == "3":
                        k2 = (un_ln(tinvdelta) * k1)
                    print("k2=" + str(k2))


            if tinvdelta == "":
                t2flag = True
                if T1 == "":
                    if T2 == "":
                        T2flag = False  # T2 was already checked
                        # tinv is being looked for
                        if arrchoice == "2":
                            tinvdelta = (Ea / R) / lnk2k1
                        if arrchoice == "3":
                            tinvdelta = lnk2k1
                        print("((1/T1)-(1/T2))=" + str(tinvdelta))

                    else:
                        # just T1 is being looked for
                        T1 = 1/(((lnk2k1*R)/Ea)+(1/T2))
                        print("T1=" + str(T1))

                if (T2 == "") and T2flag:
                    # just T2 is being looked for
                    T2 = -1/(((lnk2k1*R)/Ea)-(1/T1))
                    print("T2=" + str(T2))
        input("[enter]")


    if choice == 27:
        kc = detectpwr(input("Kc="))
        kp = detectpwr(input("Kp="))
        T = detectpwr(input("Temp="))
        n = detectpwr(input("delta_n [default=1]="))
        R = 0.08206 

        if n == "":
            n = 1

        print("*******************")
        if kc == "":
            kc = kp/((R*T)**n)
            print("Kc="+str(kc))
        if kp == "":
            kp = kc*((R*T)**n)
            print("Kp="+str(kp))

        input("[enter]")
    if choice == 28:
        kw = detectpwr(input("kw="))
        pH = detectpwr(input("ph="))
        pOH = detectpwr(input("pOH="))
        if kw == "":
            kw = pH*pOH
            print("kw="+str(kw))
        if pH == "":
            pH = kw/pOH
            print("pH="+str(pH))
        if pOH == "":
            pOH = kw/pH
            print("pOH="+str(pOH))
        input("[enter]")
    if choice == 696969696:

        # THE CODE IN THIS CHOICE IS WRONG, DO NOT USE.

        ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        leftside = []
        loopflag = True
        loopi = 0
        rxn_order = 0
        print("Enter parts in correct order...")

        # reaction order
        # ==============================================================
        # rate stays consistent                        -> zero order (power of 0)
        # moles and rate double/halve                  -> first order (power of 1)
        # moles double/halve, rate quadruples/quarters -> second order (power of 2)
        # the exponentiation of the rate from the moles determines the order

        while loopflag:
            rxn_order += rxn_order
            tmp = detectpwr(input(str(ALPHA[loopi]) + "[moles]:"))
            if tmp == "":
                loopflag = False
            else:
                leftside.append(tmp)
                loopi += 1
        print("Overall Reaction Order: " + str(1 + rxn_order))
        input("[enter]")

    # ========================================================

    if choice == 9999:
        # simplify mole ratio
        elem_mol = {}  # collect element and moles here
        highest_tenth = 1
        tmp = 0

        # kickoff
        name = detectpwr(input("enter element name [or blank to stop]: "))
        while name != "":
            mols = detectpwr(input("enter grams per mole of element      : "))

            # check how high the decimal needs to be trunc'd
            tmp_mols = mols
            nodec = float(str(tmp_mols).replace(".", ""))
            index = 0  # start at 1 to avoid div by 0
            loop_flag = True
            # there are way more efficient ways to do this but whatever
            while loop_flag:
                print(tmp_mols)
                if tmp_mols != nodec:
                    tmp_mols *= 10
                    index += 1
                if tmp_mols == nodec:
                    if index > highest_tenth:
                        highest_tenth = index
                    loop_flag = False

            # write to dict
            elem_gpm[name] = mols
            # next element (+ test for input continuation)
            name = detectpwr(input("enter element name [or blank to stop]: "))

        # once collection is done transform to a list for indexing
        elem_mol = list(elem_mol.items())

        # multiply all to the highest tenth
        for i in range(0, len(elem_mol)):
            # turn tuple into a list first
            elem_mol[i] = list(elem_mol[i])
            elem_mol[i][1] *= (10 ** highest_tenth)

        # now brute force the common factor
        loop_flag2 = True
        mod_flag = True  # true while there is modulo
        index2 = 1
        while loop_flag2:
            # check if mod flag holds up
            mod_flag = False
            for e in elem_mol:
                if (e[1] % index2) != 0:
                    mod_flag = True

                    # if modflag is not held up
            if mod_flag == False:
                # we got em
                loop_flag2 = False

        # on loop break
        print("common factor: " + str(index2))
        for e in elem_mol:
            print(str(e[0]) + str(e[1] / index2) + ":", end="")

        parts = len(elem_mol)
        input("[enter]")

    # =================================================

    if choice == 0:
        # command test
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
        # UNFINISHED
        elem_gpm = {}
        name = "PLACEHOLDER"
        gpm = 0
        while name != "":
            name = detectpwr(input("enter element name [or blank to stop]: "))
            gpm = detectpwr(input("enter grams per mole of element      : "))
            elem_gpm[name] = gpm
        parts = len(elem_gpm)
        print("")

    if choice == 42069:
        print("ok now this is epic B)")
