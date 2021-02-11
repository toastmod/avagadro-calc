from periodictable import elements

ptable = {"STP": {"Pkpa": 101.321, "Patm": 1, "Tk": 273.15, "Tf": 32, "Tc": 0, "molV": 22.4, "molV1bar": 22.7}}
for e in elements:
    ptable[e.symbol] = {"name": e.name, "mass": e.mass, "num": e.number, "ions": tuple(e.ions), "isotopes": list(e.isotopes), "charge": e.charge}
print(ptable)
print("opening file...")
dataf = open(".\\data.txt","w")
print("writing data...")
dataf.write(str(ptable))
print("closing...")
dataf.close()
print("done.\n")