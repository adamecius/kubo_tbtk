import sys
sys.path.append("/usr/local/lib/python3.6/site-packages")
import kwant
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random
welcome="Bibliotecas cargadas. ¡Bienvenidos!"

lat_c  = 0.246;
lat_vec= lat_c*np.array(((1, 0), (0.5,0.5*np.sqrt(3))));
orbs   = lat_c*np.array([(0, 0), (0, 1 / np.sqrt(3))]);
graphene = kwant.lattice.general(lat_vec, orbs);
a, b = graphene.sublattices

def create_system( L, W, sym=None, U=1.0, c=0.0, r0=(0,0) ):

  syst = kwant.Builder()
  if sym is not None:
    syst = kwant.Builder(sym)
    
  def shape(pos):
    x, y = pos;
    a0, a1 = np.linalg.inv(lat_vec).T.dot(pos);
    return ( r0[0]/lat_c <= a0 <= L/lat_c ) and ( r0[1]/lat_c <= a1 <= W/lat_c )

  def anderson(site):
      if random()*100 <= c:
        return U*(random()-0.5);
      return 0.0;

  #incorporate anderson disorder as onsites
  syst[graphene.shape(shape, r0)] = anderson;

  #incorporate hoppings
  syst[graphene.neighbors()] = -2.8;

  return syst;

def crear_cable(L, W, U=1.0, c=0.0 ):
  return create_system(L, W, sym=None, U=U, c=c );

def agregar_contactos(syst,L, W):
    tdir=-graphene.vec((1,0));
    sym = kwant.TranslationalSymmetry(tdir);
    lead = create_system(L=L, W=W, sym=sym,r0=2*tdir );
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())
    return syst;

def graficar_sistema( syst ):
  kwant.plot(syst, show=False,dpi=100);
  imp_pos = [s.pos for s,v in syst.site_value_pairs() if np.abs(v(s))>1e-5 ]
  if len(imp_pos)>0: 
    ax = plt.gca();
    ax.scatter( *np.transpose(imp_pos), c="k",zorder=3,s=1);
  plt.show();

def calcula_conductancia(fsyst,E0,nreal=1 ):
  C0 = 7.7480e-5;
  return C0*np.mean([ kwant.smatrix(fsyst, E0).transmission(1, 0) for i in range(nreal)]); 

def calcula_resistencia(fsyst,E0,nreal=1 ):

  return 1/calcula_conductancia(fsyst,E0,nreal ); 




