---
title:  "Simple is better" 
categories: 
    - Blog 
tags: 
    - Best practices 
    - PEP8 
    - Python
toc: true
---

What clean and modular code means, its importance, and some examples of best
practices.

---------------------------

## Clean and Modular Code 

The main take home message of my last two weeks working on my [GSoC_project]
can be summarized with a golden rule: write clean and modular code.


**CLEAN:** Readable, simple, and concise. <br><br>
`Clean code == Simple and clear != Redundant`
{: .notice--warning} 
{: style="text-align: center;"}


Clean code is a characteristic of production quality code that is crucial for collaboration
and maintainability in software development. Make it easy for the team to
understand and maintain.


**MODULAR:** Logically broken up into functions and modules. <br><br>
`Modular code == Broken into functions and modules =! Entangled spaghetti`
{: .notice--success}
{: style="text-align: center;"}


Modularity makes your code more organized, efficient, and reusable. This is something
that may not be evident in your first experiences working on teams, but in the
long run, it will pay off!

The importance of having clean and modular code can be summarized in three
points. First, modular code means better readability because each function will
be focused on one single purpose. Second, by having modules with a single
purpose and function, we can reuse code more easily. And third, since it is easy
to read, and to reuse, collaborating with members of your teams gets easier and
faster, speeding up development time.


Other tips that I found useful to consider while working at the [GSoC] are summarized below.

## Names of Variables

Provide descriptive names of variables. Long names are not necessarily
meaningful names. A good rule of thumb is: 


**Descriptive variables + clear variable names --> no comments** ✅ <br>
**Arbitrary variables + unclear variable names --> needs comments** ❌
{: .notice}
{: style="text-align: center;"}

-   Avoid single letter names. Exceptions can be made based on the readers of
    your code when iterables are involved, but in general try to assign the best name to variables.

    ```
    a = [15, 30, 45]        # single letter ❌

    angles = [15, 30, 45]   # descriptive ✅ 
    ```

-   Because python is untyped, specify the type of variable in the name if possible.

    ```
    is_in_list = True                 # Boolean Variable

    angle_list = [15, 30, 35]         # List
    
    dict_lipid_types = {'POPC': 15, 
                        'POPE': 10, 
                        'CHOL':5 }    # Dictionary
    ```

-   Prefer list comprehensions over explicit lists when you can create a list from an iterable.

```
angle_list = [15, 30, 45, 60, 75, 90]    # explicit list  ❌

angle_list = [i*15 for i in range(1,7)]  # list comprehension  ✅ 
```

- Concise variable names can make code easier and more effective to read.

Concise names are easy to read and names may be even intuitive or common for
other members of your team. 

For example, when working with lipid bilayers, it is very common to identify
atoms by lipid type or by atom name. An effective way to ierate over items of
`dict_lipids_ecoli` would be


```
dict_lipids_ecoli = {'PE': 34, 
                     'PG': 33, 
                     'CL': 33 }   

for lipid_type, number in dict_lipids_ecoli.items():
    print(lipid_type, number)

```


## Functions 
Each function should be focused on ONE thing. Avoid unnecessary
side effects and keep it focused on one single task.


- A _code smell_ is a function with more than 3 arguments. In many cases, using
  3 arguments may be more effective, but a good rule of thumb to keep in mind is:

**If the function has too many parameters, consider breaking into smaller functions** 
{: .notice} 
{: style="text-align: center;"}

```
function_name(param_1, ..., param_10)    # too many! ❌

function_name(param_1, ..., param_5)     # a lot, but better ⚠️

function_name(param_1, param_2, param_3) # best ✅ 
```



## Style Guide
PEP8 provides guidelines on how to write code in Python. For example:

- Use whitespace properly, and use consistent indentation (4 spaces / indent).
- Line length should be shorter than 79 characters to keeping your code human readable. 

Here an example of how select atoms using the [MDAnalysis] package following PEP8 guidelines

```
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD

u = mda.Universe(PSF, DCD)

# Line with 117 characters long ---> too long ❌
no_protein_backbone_by_residue = u.select_atoms("resname ALA or resname GLY and not (backbone or name CB)").residues

# Same selection but using a 70-character long  -->  ✅ good
protein_resids_noBB =  u.select_atoms("resname ALA GLY and not \
                                     (backbone or name CB)").residues
```


You can find more examples on code layout using PEP8 in [here]. For options to
install Autopep8, this [post]({% post_url 2021-07-12-Autopep8 %}) may help. More
useful cases for Atom selection using [MDAnalysis] are available in this
[tutorial].


In summary, the advantages of writing clean and modular code, in addition to
best practices, will take you one step closer to a more readable code, which
will be easier for your team to understand and maintain. Because simple code is
better. 


[GSoC]: https://summerofcode.withgoogle.com
[GSoC_project]: https://summerofcode.withgoogle.com/projects/#5098282306502656
[MDAnalysis]: https://github.com/MDAnalysis
[here]: https://www.python.org/dev/peps/pep-0008/?#code-lay-out
[tutorial]: https://www.mdanalysis.org/MDAnalysisTutorial/basics.html