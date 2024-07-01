
# Solbatazio energia askeen kalkulua transformazio alkimikoen bitartez

Nire Gradu Amaierako Lana osatzen duen errepositorioa da hau. Lan honetan **dinamika molekularra** eta **transformazio alkimikoak** erabiliz hainbat molekulen solbatazio energia askeak kalkulu dira. Simulazioak burutzeko beharrezkoak izan diren fitxategiak eta programak aurki daitezke hemen.

<p align="middle">
  <img src="/Molekulen_irudiak_VMD/1.png" width="100" />
  <img src="/Molekulen_irudiak_VMD/3.png" width="100" />
  <img src="/Molekulen_irudiak_VMD/7.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/11.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/21.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/23.png" width="100" />
  <img src="/Molekulen_irudiak_VMD/28.png" width="100" />
  <img src="/Molekulen_irudiak_VMD/5.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/6.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/14.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/27.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/25.png" width="100" /> 
  <img src="/Molekulen_irudiak_VMD/16.png" width="100" /> 
</p>

## Erabilera

Simulatutako 29 molekulei buruzko informazioa `Molekulak` direktorioko `Molecule.pdb` eta `Molecule.sdf` fitxategietan aurki daiteke. `Molecule.pdb` fitxategi horiek berak erabili dira `VMD` softwarearen bitartez `Molekulen_irudiak_VMD` direktorioko irudiak sortzeko.

`Programak` karpetan simulazioak egiteko programak biltzen dira. Bi zatitan egin dira simulazioak, horregatik bereizten dira `Hasierako_simulazioak` eta `Produkzio_simulazioa`. Inork programa hauek erabili nahi izanez gero, beharrezkoa izan daiteke fitxategien bide izenak aldatzea.

Lortutako emaitzen fitxategiak `Emaitzak_metrikak` direktorioko `outputs` karpetetan aurki daitezke. `npy` formatuko fitxategietan gordetzen dira, horiek baitira estatistikak kalkulatzeko `alchemical_metrikak_osoa_erroreekin.ipynb` notebookak irakurtzen dituenak. `npy` fitxategietako tentsoreetako datuen egitura ikusteko, `U_kln_lehena` eta `U_kln_azkena` `csv` fitxategiak erabil daitezke. Aipatu berri den notebook horrekin kalkulatu dira energia askeak eta metrikak. Emaitzekin sortutako irudiak `outputs/outputs_metrikak` direktorioan daude ikusgai.

## Beharrezko paketeak

Lan honetan `Python 3.11.8` erabili dira. Horrez gain, jarraian zerrendatzen dira beharrezkoak izan diren paketeak eta haien bertsioak.
```
# Izena                   Bertsioa                   
ambertools                23.3           
amberutils                21.0                    
matplotlib-base           3.8.3          
matplotlib-inline         0.1.6          
numpy                     1.26.4         
openff-amber-ff-ports     0.0.4          
openff-forcefields        2024.01.0      
openff-interchange        0.3.21         
openff-interchange-base   0.3.21         
openff-models             0.1.2          
openff-toolkit            0.15.2         
openff-toolkit-base       0.15.2         
openff-units              0.2.1          
openff-utilities          0.1.12         
openmm                    8.1.1          
openmmforcefields         0.11.2         
openmmtools               0.23.1         
pandas                    2.2.1          
scienceplots              1.0.1          
```

## Simulatutako molekulak

Simulatutako molekulen zerrenda. Hasierako simulazioetan *Errazak I* erabili dira. Bukaerakoetan, guztiak. *Errazak* eta *zailak* bereizteko, lotura birakorren kopurua kontuan hartu da (3 edo gehiago dituztenak zailtzat hartu dira). Molekula bakoitzari dagokion irudia `Molekulen_irudiak_VMD` direktorioan aurki daiteke, taulan egokitutako indizearen arabera.

<table>
    <thead>
        <tr>
            <th colspan="2">Errazak I</th>
            <th colspan="2">Errazak II</th>
            <th colspan="2">Zailak</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>1</td>
            <td>1,4-dioxanoa</td>
            <td>10</td>
            <td>Neopentanoa</td>
            <td>22</td>
            <td>Flurbiprofenoa</td>
        </tr>
        <tr>
            <td>2</td>
            <td>1-amino-4-hidroxi-9,10-Antrazenodiona</td>
            <td>11</td>
            <td>Dibentzo-p-dioxina</td>
            <td>23</td>
            <td>Ibuprofenoa</td>
        </tr>
        <tr>
            <td>3</td>
            <td>2-bromo-2-metil-propanoa</td>
            <td>12</td>
            <td>Ziklohexanoa</td>
            <td>24</td>
            <td>Ketoprofenoa</td>
        </tr>
        <tr>
            <td>4</td>
            <td>[(2S)-butan-2-il]-nitratoa</td>
            <td>13</td>
            <td>Etoxibentzenoa</td>
            <td>25</td>
            <td>2-N-etil-6-(metil-sulfonil)-4-N-(propan-2-il)-1,3,5-triazin-2,4-diamina</td>
        </tr>
        <tr>
            <td>5</td>
            <td>3-metilbutan-2-ona</td>
            <td>14</td>
            <td>Propan-1-ola</td>
            <td>26</td>
            <td>Dimetil-sulfatoa</td>
        </tr>
        <tr>
            <td>6</td>
            <td>3-fenilpropan-1-ola</td>
            <td>15</td>
            <td>1,4-dimetilpiperazina</td>
            <td>27</td>
            <td>Kloropirifosa</td>
        </tr>
        <tr>
            <td>7</td>
            <td>But-1-enoa</td>
            <td>16</td>
            <td>Metil-ziklopropano-karboxilatoa</td>
            <td>28</td>
            <td>Azido butirikoa</td>
        </tr>
        <tr>
            <td>8</td>
            <td>Ziklopentanola</td>
            <td>17</td>
            <td>Etilenoa</td>
            <td>29</td>
            <td>(2Z)-3,7-dimetilokta-2,6-dien-1-ola</td>
        </tr>
        <tr>
            <td>9</td>
            <td>Heptan-4-ona</td>
            <td>18</td>
            <td>Metanamina</td>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
        </tr>
        <tr>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
            <td>19</td>
            <td>3-kloroprop-1-enoa</td>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
        </tr>
        <tr>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
            <td>20</td>
            <td>Bromoformoa</td>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
        </tr>
        <tr>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
            <td>21</td>
            <td>Endrina</td>
            <td>&nbsp;</td>
            <td>&nbsp;</td>
        </tr>
    </tbody>
</table>

## Simulazioen emaitzak
Emaitzarik egokienak uretako sistema NPT multzoan mugalde baldintza periodikoekin eta hutsekoa NVT multzoan mugalde baldintza periodikorik gabe simulatuta lortu dira. Energiaren estimatzailerik egokiena BAR dela ikusi da. Modu honetan lortutako datuak beheko irudian alderatzen dira datu esperimentalekin.

<p align="middle">
  <img src="/Emaitzak_metrikak/outputs/outputs_metrikak/Metrikak_guztiak_BAR_2err_SCALED.png" width="600" />
</p>
<p align="middle">
    <em><b>1. irudia:</b> Lan honetan egindako simulazioen emaitzen eta balio esperimentalen arteko alderaketa. NVT-NPTnp simulazioak.</em>
</p>

Lortutako balioak [MOBLEY Lab](https://github.com/MobleyLab/FreeSolv) taldeak egindako simulazioen emaitzekin alderatu dira. Talde honek 642 molekula simulatu zituen arren, lan honetan erabilitako 29 molekulak baino ez dira irudikatu ondoko grafikoan.

<p align="middle">
  <img src="/Emaitzak_metrikak/outputs/outputs_metrikak/Metrikak_guztiak_MOBLEY.png" width="600" />
</p>
<p align="middle">
    <em><b>2. irudia:</b> MOBLEY Lab taldeak molekula berak simulatzean lortutako balioak.</em>
</p>

1. irudia sortzeko erabilitako datuak, esperimentalak zein simulazioen emaitzak, `emaitzak_metrikak` direktorioko [`emaitzen_laburpena_NPT_NVTnp.csv`](Emaitzak_metrikak/emaitzen_laburpena_NPT_NVTnp.csv) fitxategian aurki daitezke.
   
## Lizentzia

```CC BY 4.0```
