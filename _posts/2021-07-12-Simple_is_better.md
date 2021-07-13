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

The main take home message of my last two weeks working on my [project]
can be summarized with a golden rule: write clean and modular code.


**CLEAN:** Readable, simple, and concise. <br><br>
`Clean code == Simple and clear != Redundant`
{: .notice--warning} 
{: style="text-align: center;"}


Clean code: a characteristic of production quality code that is crucial for collaboration
and maintainability in software development. Make it easy for the team to
understand and maintain.


**MODULAR:** Logically broken up into functions and modules. <br><br>
`Modular code == Broken into functions and modules =! Entangled spaghetti`
{: .notice--success}
{: style="text-align: center;"}


Modularity: It makes your code more organized, efficient, and reusable. This is something
that may not be evident in your first experiences working on teams, but in the
long run, it will pay off!

The importance of having clean and modular code can be summarized in three
points. First, modular code means better readability of the code. Second, by
having modules, we can reuse code more easily. And third, since it is easy to
read, and to reuse, collaborating with members of your teams get easier and
faster, speeding up development time.

However, more modules in your code is not always more effective. There is always
a balance. 

Imagine we have a trajectory and a coordinates file that we obtained from
Molecular Dynamics (MD) simulations. All we know is that the file contains a
complex membrane, and we need to find how many lipid types we have in the
system, the number of atoms (`n_atoms`) and. The `MDAnalysis` package is
probably the best option to start exploring the system. 



Other tips that I found useful to consider while workin on GSoC are summarized below:

## Names of Variables

Provide descriptive names of variables. Long names are not necessarily
meaningful names. A good rule of thumb is: 


**Explicit variables + clear variable names --> no comments** ✅ <br>
**Implicit variables + unclear variable names --> needs comments** ❌
{: .notice}
{: style="text-align: center;"}

-   Avoid single letter names.  (Exceptions can be made based the readers of
    your code)

    ```
    a = [15, 30, 45]        # single letter ❌

    angles = [15, 30, 45]   # descriptive ✅ 
    ```

-   Because python is untyped Specify type of variable in the name if possible.

    ```
    is_in_list = True                 # Boolean Variable

    angle_list = [15, 30, 35]         # List
    
    dict_lipid_types = {'POPC': 15, 
                        'POPE': 10, 
                        'CHOL':5 }    # Dictionary
    ```

-   Prefer list comprehensions over explicit lists. 

```
angle_list = [15, 30, 45, 60, 75, 90]    # explicit list  ❌

angle_list = [i*15 for i in range(1,6)]  # list comprehension  ✅ 
```

- Concise variable names can be more effective in certain functions.

```
for lipid_type, number in dict_lipid_types.items():
    print(lipid_type, number)
```


## Functions 
Each function should be focused on ONE thing. Avoid unnecessary
side effects and keep it focused on one single task.




- A _code smell_ for functions with more than 3 arguments. Although
doesn’t apply for 100% of the cases, using 3 arguments may be more
effective.  A rule of thumb to keep in mind is, 

**If the function has too many parameters, <br> consider to break the function into more modules.**
{: .notice}
{: style="text-align: center;"}

```
function_name(param_1, ..., param_10)    # not very good ❌

function_name(param_1, ..., param_5)     # better ⚠️

function_name(param_1, param_2, param_3) # best ✅ 
```



## Style Guide
PEP8 provides guidelines on how to write code in Python. For example,
- Use whitespace properly and use consistent indentation (4 spaces / indent).
- Line length should be shorter than 79 characters. 


```
lipid_list_ecoli = ['CL', 'PE', 'PG' ]
common_lipid_types = ['PC', 'PE', 'PS']
not_very_common_lipid_types = ['PI', 'TAG', 'PUFA]

for lipid in lipid_list_ecoli:
    if lipid is in dict_lipid_types.keys() or is in dict_not_very_common_lipid_types.keys():
        print( lipid )
        
```
    
You can find more examples on code lay-out using PEP8 in [here]. For options
to install Autopep8, this [post]({% post_url 2021-07-12-Autopep8 %}) may
help.


Following the practices here mentioned will take you one step closer to a more
readable code, which will be easier for your team to understand and maintain. 

In summary, write clean and modular code. Because simpler is better. 


[GSoC]: https://summerofcode.withgoogle.com
[project]: https://summerofcode.withgoogle.com/projects/#5098282306502656
[MDAnalysis]: https://github.com/MDAnalysis
[here]: https://www.python.org/dev/peps/pep-0008/?#code-lay-out
