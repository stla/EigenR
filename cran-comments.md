## Release summary

This is a resubmission. The first submission has been archived because of a 
R CMD CHECK warning on Debian:

* checking whether package ‘EigenR’ can be installed ... [297s/297s] WARNING
Found the following significant warnings:
  EigenR.cpp:10:44: warning: imaginary constants are a GCC extension

It is due to the presence of `1i` in this code, which is the "imaginary
constant in question:

```
return Re.cast<std::complex<double>>() + 1i * Im.cast<std::complex<double>>();
```

So I replaced this code with:

```
const std::complex<double> I_ {0.0, 1.0};
return Re.cast<std::complex<double>>() + I_ * Im.cast<std::complex<double>>();
```

I tried to check on r-hub but I can't because of issues with my internet 
connection.


## Test environments

* ubuntu 18.04, R 3.6.3
* win-builder (devel)


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
