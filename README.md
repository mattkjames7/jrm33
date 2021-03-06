# jrm33

JRM33 model (Connerney et al. 2022) implementation using Python.

## Installation

Install using `pip`:

```bash
pip3 install jrm33 --user
```

Or by cloning this repo:

```bash
git clone https://github.com/mattkjames7/jrm33.git
cd jrm09

#EITHER create a wheel and install (replace X.X.X with the version number):
python3 setup.py bdist_wheel
pip3 install dist/jrm33-X.X.X-py3-none-any.whl --user

#OR install directly using setup.py
python3 setup.py install --user
```

## Usage

The model accepts right-handed System III coordinates either in Cartesian form (`jrm33.ModelCart()`) or in spherical polar form (`jrm33.Model()`), e.g.:

```python
import jrm33

#get some Cartesian field vectors (Deg keyword is optional)
Bx,By,Bz = jrm33.ModelCart(x,y,z,Deg=13)

#or spherical polar ones
Br,Bt,Bp = jrm33.Model(r,theta,phi,Deg=13)
```

Please read the docstrings for `jrm33.Model()` and `jrm33.ModelCart()` using `help` or `?` e.g. `help(jrm33.Model)` .

There is also a test function which requires `matplotlib` to be installed:

```python
#evaluate the model at some R
jrm33.Test(R=0.85)
```

which produces this (based on figure 4 of Connerney et al. 2018):

![jrm09test.png](jrm33test.png)

## References

Connerney, J. E. P., Kotsiaros, S., Oliversen, R. J., Espley, J. R.,  Joergensen, J. L., Joergensen, P. S., et al. (2018). A new model of Jupiter's magnetic field from Juno's first nine orbits. Geophysical Research Letters, 45, 2590– 2596. https://doi.org/10.1002/2018GL077312

Connerney, J. E. P., Timmins, S., Oliversen, R. J., Espley, J. R., Joergensen, J. L., Kotsiaros, S., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's Prime Mission. *Journal of Geophysical Research: Planets*, 127, e2021JE007055. [A New Model of Jupiter's Magnetic Field at the Completion of Juno's Prime Mission - Connerney - 2022 - Journal of Geophysical Research: Planets - Wiley Online Library](https://doi.org/10.1029/2021JE007055)
