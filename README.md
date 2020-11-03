# avagadro-calc
A python script i wrote to breeze through my university chemistry course.

# installation
Please install [periodictable](https://github.com/pkienzle/periodictable) before use.
`pip install periodictable`
**This will require you to use python 3.7 or earlier** as 3.8 dropped support for this library.
If you want to use the script with 3.8 you will have to delete the lines up to the comment that says `# CUT HERE`
This means you cant use periodic table features.

# usage
While the interface is pretty straight forward, there are a few hidden tricks.
## using *10^
When entering a number, you can end it with "e" to prompt a "*10^" entry.
This will (obviously) multiply the number by 10 to number entered.
example:
```
[#]: 3
===================================
atoms: 7e
*10^6
moles: 1.1624045167718366e-17
```
## using getting periodic table data
Currently only mass and atomic number are available, but if you know how to use the periodictable you can add the other data to the dictionary generator in the header of the script.
You can enter a mass number like this (in this case Helium):
```
He:mass
```
or an atomic number like this:
```
He:num
```
in the 0 command, you can output the name of an element:
```
[$]>>He:name
helium
[enter]
```



