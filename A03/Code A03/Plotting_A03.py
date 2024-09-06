#!/usr/bin/env python
# coding: utf-8

# ## Plotting Graphs for A03
# ### Nirmal S. [234103107]

# In[30]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[31]:


temp_exact = pd.read_csv('5. Temp (Analytical).csv', header=None).iloc[:,:-1].iloc[::-1].reset_index(drop=True)
temp_cg = pd.read_csv('1.2 Temp (CG).csv', header=None).iloc[:,:-1].iloc[::-1].reset_index(drop=True)
temp_pcg_J = pd.read_csv('2.2 Temp (PCG_Jacobi).csv', header=None).iloc[:,:-1].iloc[::-1].reset_index(drop=True)


# In[ ]:


## Function to plot the contour plot

def contour_plot(df, text):
    
    font1 = {'family':'Times New Roman','color':'Black','size':14}
    
    x = np.arange(len(df.columns))
    y = np.arange(len(df.index))
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots()
    cs = ax.contourf(X, Y, df.values, cmap='coolwarm')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    
    plt.title("Temperature Contour "+text, fontdict=font1)
    
    plt.colorbar(cs)
    plt.savefig("Temperature Contour "+text+".jpg")
    plt.show()


# In[ ]:


## Function to find the series of at midsection

def horizontal_midsection(df):
    if ((len(df.columns)-1)%2==0):
        y = df[int((len(df.columns)-1)/2)]
    else:
        y = 0.5*(df[int(len(df.columns)/2)] + df[int(len(df.columns)/2-1)])

    return y

def vertical_midsection(df):
    if ((len(df.index)-1)%2==0):
        x = df.loc[int((len(df.columns)-1)/2)]
    else:
        x = 0.5*(df.loc[int(len(df.columns)/2)] + df.loc[int(len(df.columns)/2-1)])

    return x


# In[ ]:


## Function to plot temp values along midsections

def plot_horizontal_mid(df, text, i):
    plt.figure(figsize=(6, 6))

    font1 = {'family':'Times New Roman','color':'Black','size':20}
    font2 = {'family':'Times New Roman','color':'Black','size':14}

    colours = ['black', 'red', 'orange', 'yellow', 'green', 'blue', 'violet']

    y = np.linspace(0, 1, len(temp_exact.columns))
    x_exact = horizontal_midsection(temp_exact)
    x_df = horizontal_midsection(df)

    step = 5
    plt.plot(x_exact[::step], y[::step], label='Temp_analytic', marker='^', color=colours[0])
    plt.plot(x_df, y, label=text, color=colours[i])

    plt.xlim(0,1)
    plt.ylim(0,1)

    plt.legend()

    plt.title("Temperation along y-axis at 0.5 x", fontdict=font1)

    plt.xlabel("Temperature at x=0.5", fontdict=font2)
    plt.ylabel("Y", fontdict=font2)
    
    plt.savefig("Horizontal_Mid_"+text+".jpg")
    plt.show()
    
    
def plot_vertical_mid(df, text, i):
    plt.figure(figsize=(6, 6))

    font1 = {'family':'Times New Roman','color':'Black','size':20}
    font2 = {'family':'Times New Roman','color':'Black','size':14}

    colours = ['black', 'red', 'orange', 'yellow', 'green', 'blue', 'violet']

    x = np.linspace(0, 1, len(temp_exact.columns))
    y_exact = vertical_midsection(temp_exact)
    y_df = vertical_midsection(df)

    step = 5
    plt.plot(x[::step], y_exact[::step], label='Temp_analytic', marker='^', color=colours[0])
    plt.plot(x, y_df, label=text, color=colours[i])

    plt.xlim(0,1)
    plt.ylim(0,0.3)

    plt.legend()

    plt.title("Temperation along x-axis at 0.5 y", fontdict=font1)

    plt.xlabel("X", fontdict=font2)
    plt.ylabel("Temperature at y=0.5", fontdict=font2)
    plt.savefig("Vertical_Mid_"+text+".jpg")
    plt.show()


# In[ ]:


## Temperature Contour

contour_plot(temp_exact, 'Analytic')
contour_plot(temp_cg, 'Conjugate Gradient')
contour_plot(temp_pcg_J, 'Preconditioned CG - Jacobi')


# In[ ]:


## Temperature along 'y' about horizontal midsection

plot_horizontal_mid(temp_cg, 'Temp_Conjugate-Gradient', 1)
plot_horizontal_mid(temp_pcg_J, 'PreConditioned-CG-Jacobi', 5)


# In[ ]:


## Temperature along 'x' about vertical midsection

plot_vertical_mid(temp_cg, 'Temp_Conjugate-Gradient', 1)
plot_vertical_mid(temp_cg, 'PreConditioned-CG-Jacobi', 5)


# In[ ]:





# In[ ]:




