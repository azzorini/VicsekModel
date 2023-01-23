# VicsekModel

In this repository I give some code to simulate the Vicsek model. All the code is written in C++ except for a Jupyter Notebook that was used to generate GIFs from the data.

The list of programs included in this repository is the following:

1. **vicsek_gif.cpp:** Generates the data to produce a GIF of the Vicsek model with the given parameters.
2. **vicsek_gif_linklist.cpp:** Same as the previous one but it uses a link list to optimize the code.
3. **BirdsGif.ipynb:** Jupyter notebook to produce the GIFs.
4. **vicsek_va.cpp:** Generates a plot for the order parameter as a function of time. It uses a link list to optimize the code.
5. **vicsek_va_tra.cpp:** Generates a plot for the order parameter as a function of the noise. It uses link lists to optimize the code.
6. **vicsek_va_without_link.cpp:** Same thing than the previous one but it does not use link lists. This makes the code simpler but also slower.
