import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter

# PLOT 1 sigmaT vs Q2 avec W fixé
x, y = np.loadtxt("/Users/mr282803/test/cs_test.txt", unpack=True)

plt.plot(x, y, color='red', linestyle='--', label='W = 2.5 GeV')
plt.xscale('log')
plt.yscale('log')

ax = plt.gca()

class LogFormatterPlain(LogFormatter):
    def __call__(self, x, pos=None):
        
        if x in [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]:
            return f"{x:g}"
        else:
            return ""  

ax.xaxis.set_major_formatter(LogFormatterPlain())
ax.yaxis.set_major_formatter(LogFormatterPlain())

ax.tick_params(direction='in', width=2, length=6, which='major')
ax.tick_params(direction='in', width=1.5, length=3, which='minor')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

ax.set_xlabel(r"$\mathbf{Q^{2}\ [GeV^{2}]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{\sigma_T\ [nb]}$", fontsize=14)
ax.set_title(r"$\mathbf{Transverse\ \sigma_T\ (\gamma^{*} + p \rightarrow \phi + p)}$", fontsize=15)
ax.legend(fontsize=12)

ax.set_xlim(0.1, 100)
ax.set_ylim(0.01, 1000)

plt.tight_layout()
plt.show()








# PLOT 2 sigmaT vs W avec Q2 fixé

x, y = np.loadtxt("/Users/mr282803/test/cs_test2.txt", unpack=True)

plt.plot(x, y, color='red', linestyle='--', label=r"$Q^{2} = 2.5\ \mathrm{GeV}$")
plt.xscale('log')
plt.yscale('log')

ax = plt.gca()

class LogFormatterPlain(LogFormatter):
    def __call__(self, x, pos=None):
        
        if x in [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]:
            return f"{x:g}"
        else:
            return ""  

ax.xaxis.set_major_formatter(LogFormatterPlain())
ax.yaxis.set_major_formatter(LogFormatterPlain())

ax.tick_params(direction='in', width=2, length=6, which='major')
ax.tick_params(direction='in', width=1.5, length=3, which='minor')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

ax.set_xlabel(r"$\mathbf{W\ [GeV]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{\sigma_T\ [nb]}$", fontsize=14)
ax.set_title(r"$\mathbf{Transverse\ \sigma_T\ (\gamma^{*} + p \rightarrow \phi + p)}$", fontsize=15)
ax.legend(fontsize=12)

ax.set_xlim(1, 100)
ax.set_ylim(1, 100)

plt.tight_layout()
plt.show()





# PLOT 3 R vs Q2

x, y = np.loadtxt("/Users/mr282803/test/cs_test3.txt", unpack=True)

plt.plot(x, y, color='red', linestyle='--')
plt.xscale('log')
plt.yscale('log')

ax = plt.gca()

class LogFormatterPlain(LogFormatter):
    def __call__(self, x, pos=None):
        
        if x in [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]:
            return f"{x:g}"
        else:
            return ""  

ax.xaxis.set_major_formatter(LogFormatterPlain())
ax.yaxis.set_major_formatter(LogFormatterPlain())

ax.tick_params(direction='in', width=2, length=6, which='major')
ax.tick_params(direction='in', width=1.5, length=3, which='minor')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

ax.set_xlabel(r"$\mathbf{Q^{2}\ [GeV^{2}]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{R = \frac{\sigma_{L}}{\sigma_{T}}\ [nb]}$", fontsize=14)
ax.set_title("Ratio R", fontsize=15)
ax.legend(fontsize=12)

ax.set_xlim(0.1, 100)
ax.set_ylim(0.01, 10)

plt.tight_layout()
plt.show()





# PLOT 4 dsimgaTdt + episolon*dsigmaLdt  vs t avec Q2 et W fixé

x, y = np.loadtxt("/Users/mr282803/test/cs_test4.txt", unpack=True)

plt.plot(x, y, color='red', linestyle='--', label=r"$Q^{2} = 2.5\ \mathrm{GeV^{2}}\ \text{and}\ W = 2.5\ \mathrm{GeV}$")

plt.yscale('log')

ax = plt.gca()

class LogFormatterPlain(LogFormatter):
    def __call__(self, x, pos=None):
        
        if x in [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]:
            return f"{x:g}"
        else:
            return ""  

ax.yaxis.set_major_formatter(LogFormatterPlain())

ax.tick_params(direction='in', width=2, length=6, which='major')
ax.tick_params(direction='in', width=1.5, length=3, which='minor')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

ax.set_xlabel(r"$\mathbf{t_{min} - t\ [GeV^{2}]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{\frac{d\sigma}{dt}\ (\gamma^{*} + p \rightarrow \phi + p)\ [\mathrm{nb/GeV^2}]}$", fontsize=14)

ax.legend(fontsize=12)

ax.set_xlim(0.0, 2.0)
ax.set_ylim(0.1, 100)

plt.tight_layout()
plt.show()





# PLOT 5 cross section total dsigma / dQ^2 dt xb partie emission du photon QED * gamma* + p -> phi + p (QCD) ce qui done bien ep->e'p'phi avec xb (W) et Q2 fixé

x, y = np.loadtxt("/Users/mr282803/test/cs_test5.txt", unpack=True)

plt.plot(x, y, color='red', linestyle='--', label=r"$Q^{2} = 2.5\ \mathrm{GeV^{2}}\ \text{and}\ W = 2.5\ \mathrm{GeV}$")

plt.yscale('log')

ax = plt.gca()

class LogFormatterPlain(LogFormatter):
    def __call__(self, x, pos=None):
        
        if x in [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]:
            return f"{x:g}"
        else:
            return ""  

ax.yaxis.set_major_formatter(LogFormatterPlain())

ax.tick_params(direction='in', width=2, length=6, which='major')
ax.tick_params(direction='in', width=1.5, length=3, which='minor')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

ax.set_xlabel(r"$\mathbf{t_{min} - t\ [GeV^{2}]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{\frac{d\sigma}{dtdxbdQ^2}\ (e + p \rightarrow \phi + p)\ [\mathrm{nb/GeV^4}]}$", fontsize=14)


ax.legend(fontsize=12)

ax.set_xlim(0.0, 2.0)
ax.set_ylim(0.0001, 10)

plt.tight_layout()
plt.show()



# PLOT 6 TAU

x, y = np.loadtxt("/Users/mr282803/test/cs_test6.txt", unpack=True)

plt.plot(x, y, color='red', linestyle='--', label=r"$Q^{2} = 2.5\ \mathrm{GeV^{2}}\ \text{and}\ W = 2.5\ \mathrm{GeV}$")

plt.yscale('log')

ax = plt.gca()

class LogFormatterPlain(LogFormatter):
    def __call__(self, x, pos=None):
        
        if x in [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]:
            return f"{x:g}"
        else:
            return ""  

ax.yaxis.set_major_formatter(LogFormatterPlain())

ax.tick_params(direction='in', width=2, length=6, which='major')
ax.tick_params(direction='in', width=1.5, length=3, which='minor')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

ax.set_xlabel(r"$\mathbf{t_{mTHETATHETATHETAin} - t\ [GeV^{2}]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{\frac{d\sigma}{dtdxbdQ^2}\ (e + p \rightarrow \phi + p)\ [\mathrm{nb/GeV^4}]}$", fontsize=14)


ax.legend(fontsize=12)

ax.set_xlim(0.0, 200)
ax.set_ylim(0.0001, 1000)

plt.tight_layout()
plt.show()


