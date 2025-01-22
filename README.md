# Revisiting RSA-polynomial problem and semiprime factorization

## Introduction

This is a Python implementation of lattice-based factoring attacks proposed in **Revisiting RSA-polynomial problem and semiprime factorization**[^RSAPoly].

## Requirements

- [**SageMath**](https://www.sagemath.org/) 9.5 with Python 3.10

You can check your SageMath Python version using the following command:

```commandline
$ sage -python --version
Python 3.10.12
```

Note: If your SageMath Python version is older than 3.9.0, some features in given scripts might not work.

## Usage

The standard way to run the attack with the specific parameters requires passing them as command line arguments `sage -python attack.py <prime_bit_length> <x_bit_length> <m> <theorem: BM or JM>`. For instance, to run the attack with $\ell=256$, $\xi=85$, $m=4$ and *theorem BM*, please run `sage -python attack.py 256 85 4 BM`:

```commandline
RSAPoly$ sage -python attack.py 256 85 4 BM
The parameters:
e0 = -2641526499254901730039349217926516765616381889365300726452178305804358863736
e1 = 57438931253243610422613711846019068095050266700370703667866622032019712740873
e2 = 66266325610106000626854456229139997576312103656593797866926097881527586607327
Found primes:
p = 66266325610106000626854456229139997576312103656594016719405703683267927490583
q = 63624799110851098896815107011213480810695721767228287011952055263786035407527
The attack costs 1.622 seconds...
```

For instance, to run the attack with $\ell=512$, $\xi=175$, $m=4$ and *theorem JM*, please run `sage -python attack.py 512 175 4 JM`:

```commandline
RSAPoly$ sage -python attack.py 512 175 4 JM
The parameters:
e0 = -5068226797795168185036327209560464514961121601551008519753233044488905226664913792710770806576448450304195219643464878599517226385728134594560578871377836
e1 = 6388417077923910472726193793241470703715312574608942669651980058637483919581399033322445295752891263031619007535980771479468450280554966278255748859800913
e2 = 12195216458738193520114447207311002639072207655438355402226099956945039637941252861336657785044222624017923400479301292130516334362443784910677442429983747
Found primes:
p = 12195216458738193520114447207311002639072207655438355402226099956945039637941252861336657785044222624395961322102450238754362907109216973312437808620587419
q = 7126989660943025335078119997750538124111086053887346882472866912456134411276339068625886978467774173492799559522021581669875624584080669742385910440979959
The attack costs 4.718 seconds...
```

## Notes

All the details of the numerical attack experiments are recorded in the `attack.log` file.

[^RSAPoly]: Zheng M., "Revisiting RSA-polynomial problem and semiprime factorization" | [PDF](https://mengcezheng.github.io/docs/Zheng24b.pdf)
