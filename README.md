
# StaticDESim

Simulation and optimization of distributed energy system under rated condition.

![](./img/1.jpg)

`StaticDESim` provides simple models for some thermal engines and heat
transfers under rated condition. It can help design and optimize distributed
energy system, which can provide electricity, heating and cooling at the same
time.

The object is a grid-connected distributed energy system, called Green Energy
Island (GEI), composed of a variety of energy technologies. It satisfies the
need of people and equipments in a remote island by supplying electricity, hot
steam, and cooling. 

- There are three facilities, Micro Gas Turbine (MGT), Aqueous Lithium-Bromine
  Single-Effect Absorption Chiller (AC_ALB), and R123 Organic Recycle Cycle
  (ORC_R123), providing the majority of three products. Nature gas or other
  substituted fuel burns in MGT to produce electricity, then some waste heat
  is recovered by AC_ALB and ORC_R134a. One produces cool water, and the other
  one produces electricity.
- If the load of hot steam and cooling cannot be met, Heat Pump (HP) and R134a
  Vapour Compression Chiller (VCC_R134a) can be powered to supply hot steam and
  cooling, respectively. More electricity can be bought from the external grid.

---

The code resulted from my bachelor thesis. I did not know too much about coding
that time, and it has not been under active development since then. Hope it can
provide a starting point somehow.
