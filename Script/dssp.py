""" Script de l'implementation de la methode DSSP
    pour l'assignation des structures secondaires
    des proteines
----------------------------------------------------------------
BARRY Fatoumata Binta
AHMAR-ERRAS Noura
M2 BI Université Paris Diderot
Année : 2020-2021
----------------------------------------------------------------

"""
# Modules utiles
import sys
import os
import math

# Fonctions


def extraction_atome(fichier_pdb):
    """ Lit un fichier pdb et renvoie les coordonées de tous les atomes
        impliqués dans la liaison peptidique (C, O, H, N).

    Parametres
    ----------
    fichier_pdb: fichier pdb contenant les coordonnées des atomes
    
    Returns
    -------
    dico: dictionnaire avec comme clé la position du residu dans la proteine
         et comme valeur une liste de tuples:
            - res: code 3 lettres du residus (ALA par exemple)
            - at: atome (C, H, O ou N)
            - x, y, z: coordonnees x, y, z de l'atome en Angstrom

    """
    nom = ["C","O","N","H"]
    with open(fichier_pdb, "r") as pdb:
        dico = {}
        for ligne in pdb:
            if ligne.startswith("ATOM") and ligne[12:16].strip() in nom and (ligne[21:23].strip()=="A"or ligne[21:23].strip()==""):
                posi = ligne[22:27].strip()
                at = ligne[12:16].strip()
                res = ligne[17:20].strip()
                x = ligne[32:38]
                y = ligne[40:46]
                z = ligne[48:54]
                if posi in dico:
                    dico[posi].append((res, at, x, y, z))
                else:
                    dico[posi] = [(res, at, x, y, z)]
    return dico
               

def distance(res1, res2):
    """ Permet de calculer la distance entre deux points
        ayant leurs cooronnées x, y et z

    Parametres
    ----------
    res1, res2: tuple de type (x, y, z)
    
    Returns
    -------
    distance: float, retourne la distance entre les deux points

    """
    x1, y1, z1 = float(res1[0]), float(res1[1]), float(res1[2])
    x2, y2, z2 = float(res2[0]), float(res2[1]), float(res2[2])
    norme = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    distance = math.sqrt(norme)
    return distance


def liaison_hydrogene(dico):
    """ Permet de determiner l'existence d'une liaison
       hydrogene (LH) entre le CO du residu i et le NH du i+n,
       ou entre le NH de i et le CO de i-n, "n" allant de 1 au 
       nombre total de residu de la chaine,en se basant sur 
       l'energie electrostatistique (E) et en appliquant la 
       condition: presence de LH si E < 0.5kcal/mol

    Parametres
    ----------
    dico: dictionnaire issue de la fonction extraction atome
        par exemple dico[1]=[("ALA", "C", "0.1", "-0.5", "22")]
    
    Returns
    -------
    liaison: dico avec comme clé la position de la paire de residus
        impliqués dans la LH et comme valeur l'energie (<0.5)
        liaison[("1,2")]=-4.99 par exemple

    """
    cle = list(dico.keys())
    val = list(dico.values())
    # q1, q2 et f ont ete donnés dans l'article, de meme que
    # la formule de l'energie
    q1 = 0.42 # charge partielle
    q2 = 0.28 # charge partielle
    f = 332 # facteur dimensionnel
    liaison = {}
    for i in range(len(cle) - 1):
        res1 = dico[cle[i]]
        C = res1[1][2:] # coordonnees de C
        O = res1[2][2:]
        
        for j in range(len(cle)):
            res2 = dico[cle[j]]
            N = res2[0][2:]
            if len(res2) == 4: # permet d'eliminer les aa sans hydrogene
                H = res2[3][2:]
                ron = distance(O,N) # distance entre l'atome "O" et "N"
                rch = distance(C,H)
                roh = distance(O,H)
                rcn = distance(C,N)
                energie = round(q1*q2*f*(1/ron + 1/rch - 1/roh - 1/rcn), 2)
                if energie <-0.5: # condition pour une LH
                    liaison[(cle[i],cle[j])] = energie
    return liaison


def pattern(liaison_hydroge):
    """ Permet d'attribuer les differents motifs, selon
        la position des atomes impliqués dans la LH
        3<=|i-i+n|<=5: il s'agit de turns
        |i-i+n|<3: pas de stucture, d'apres notre comprehension
        de l'article ("None")
        |i-i+n|>5: liaison eloignée, suspicion de brin beta
    Parametres
    ----------
    liaison hydroge: dictionnaire issue de la fonction liaison_hydrogene
    
    Returns
    -------
    pattern: dico avec comme clé le type de pattern et comme valeur
            une liste de tuples des couples de residus ayant ce pattern
            par exemple 4turns': [(3, 7), (4, 8)], les residus 3 et 4
            ont des turns de type 4turns, de meme que 4,8..

    """
    pat = {}
    for i in liaison_hydroge:
        val = liaison_hydroge[i]
        
        res1 = int(i[0])# on transforme en int pour faciliter la suite
        res2 = int(i[1])
        diff = abs(res1-res2) # calcul de la difference entre les positions
        if diff > 5:
            typ = "brin"
        elif 3 <= diff <= 5:
            typ = str(diff)+"turns"
        else:
            typ = "None"
        
        # mise en place du dico pattern
        if typ in pat:
            pat[typ].append((res1, res2))
        else:
            pat[typ] = [(res1, res2)]
    return pat


def structure(pat):
    """ Permet d'assigner le type de stucture secondaire
    helice (H) ou feuillet beta (B), suivant la localisation
    de la lh et son nombre de repetition (succession) 
    Parametres
    ----------
    pat: dictionnaire issue de la fonction pattern
    
    Returns
    -------
    struc: dico avec comme clé la position de la paire de residus
        impliqués dans la LH et comme valeur H, B ou None

    """
    struc = {}
    for i in pat:
        val = pat[i] # on recupere la liste de tuples pour chaque clé
        for j in range(len(val) - 1):
            k = val[j] # couple i
            l = val[j+1] # couple i+1 , pour voir les successions
            if "turn" in i: # pour generaliser tout en helice
                if k[0] == l[0] - 1: # on verifie qu'il y a une succession entre les deux couples
                    struc[k] = "H"
                    struc[l] = "H"
                elif (k[0] != l[0] - 1) and l not in struc:
                        struc[l] = i
            elif i == "brin": # pour statuer sur l'existance des feuillets
                if l[0] - 2 <= k[0] <= l[0] + 2: # On prend un intervalle de 2, pour de resultats coherents avec l'article
                    struc[k] = "B" # feuillet beta (B)
                    struc[l] = "B"
                elif l[0] - 2 >= k[0] >= l[0] + 2 and l not in struc:
                    struc[l] = i
            else: # pour traiter le none
                struc[k] = i
                if l not in struc:
                    struc[l] = i
    return struc


def fichier_dssp(recap, struc):
    """ Permet de generer un fichier de sortie avec 2 colonnes:
        le couple de residu, et le type de structure secondaire
        dans lequel il est impliqué (ici en l'occurence H et B et nturns) 
    Parametres
    ----------
    recap: nom du fichier de sortie
    struc: dico issu de la fonction structure
    
    Returns
    -------
    recap

    """
    with open(recap,"w") as rec:
        for i in struc:
            if struc[i] != "None": # on elimine tous les couples dont la structure n'est pas definie
                structure = struc[i] # on recupere le type de structure secondaire (H ou B ou nturn)
                rec.write("{}\t{}\n".format(i, structure)) # on ecrit dans un fichier 
    return recap


# Programme principal

if __name__ == '__main__':

    if len(sys.argv) != 3:
        sys.exit("ERREUR : il faut deux autres arguments:\n un fichier.pdb issu de reduce et un nom de fichier de sortie comme fichier.txt par exemple.")

    file_pdb = sys.argv[1]
    file_text = sys.argv[2]

    if not os.path.exists(file_pdb):
        sys.exit("ERREUR: {} n'existe pas dans ce repertoir".format(file_pdb))

    # Declaration des variables 
     
    dico_coordonnee = extraction_atome(file_pdb)
    liaison = liaison_hydrogene(dico_coordonnee)
    patt = pattern(liaison)
    struct_second = structure(patt)
    fichier_dssp(file_text, struct_second)
    
                    
           
               
                

























