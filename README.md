# FinancePy
FinancePy is a library of functions for valuing, selecting and risk-managing financial assets. 

The aim of this library for me has been to provide a comprehensive and accessible Python library for financial calculations that can be used by my students at EDHEC to learn about financial derivatives. I believe it can also be used by academics and practitioners to perform the pricing and risk-management of complex financial products. 

I intend that subsequent versions will also include asset selection, portfolio-level risk management, regulatory calculations and market analysis tools. 

In general my objectives have been:
1) To make the code as simple as possible so that students and those with a basic Python fluency can understand and check the code.
2) To keep all the code in Python so users can look through the code to the lowest level. 
3) To offset the performance impact of (2) by leveraging Numba to make the code as fast as possible without resorting to Cython.
4) To make the design product-based rather than model-based so someone wanting to price a specific exotic option can easily find that without having to worry too much about the model – just use the default – unless they want to. 
5) To make the library as complete as possible so a user can find all their required finance-related functionality in one place. This is better for the user as they only have to learn one interface.
6) To avoid complex designs as I do not want to make it too hard for unskilled Python programmers to use the library.
7) To have good documentation and easy-to-follow examples.
8) To make it easy for interested parties to contribute.

In many cases the valuations should be close to if not identical to those produced by financial systems such as Bloomberg. However for some products, larger value differences may arise due to differences in date generation and interpolation schemes. Over time I expect to reduce the size of such differences.

### Dependencies
FinancePy depends on Numpy and Numba and Scipy.

### Installation
FinancePy can be installed via pip as 
>>>> pip install financepy

### Changelog
See the changelog for a detailed history of changes 

### Contributions
Contributions are welcome, as long as you don't mind camel case ;-)

### License
MIT
