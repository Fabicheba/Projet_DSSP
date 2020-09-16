# **Implementation de l'algorithme de DSSP: assignation de structure secondaire des proteines**
### le script dssp.py permet d'assigner les structures secondaires à partir d'un fichier pdb. Pour l'utiliser, il faudra suivre les étapes suivantes. nous utiliserons le fichier pdb 3pti.pdb pour illustrer nos commandes
#### le projet DSSP est composé de 5 repertoires (Bin, Documents, Données, Script, Resultat). Se placer dans le repertoire Projet_DSSP pour executer les commandes ci-dessous
- Charger l'environ l'enviroment conda avec la commande 

``` conda env create -f Document/dssp_conda.yml```
- Charger l'environ l'enviroment conda avec la commande 

- Télécharger le fichier pdb sur le site https://www.rcsb.org/structure/3pti ou le recuperer directement dans le repertoire Donnees/3pti.pdb
-  Passer le fichier dans le programme reduce pour ajouter les hydrogenes avec la commande
``` reduce -build Donnees/3pti.pdb Donnees/3pti_h.pdb```
- Lancer le script dssp.py avec la commande
```python Script dssp.py Donnees/3pti_h.pdb Resultats/3pti_predit```
le 3pti_predit etant le fichier de sortie contenant les positions des couples de residus participant à la mise en place de la structure secondaire assignée
Nos resultats ont stockés dans le repertoire Resultats.
- Nous avons egalement utilisé le programme de dssp pour pouvoir evaluer nos resultats. Le programme est installé dans l'environnement conda, pour l'utiliser il suffit de taper la commande
```mkdssp -i Donnees/3pti_h.pdb -o Resultats/3pti_dssp ```