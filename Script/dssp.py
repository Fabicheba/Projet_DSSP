import sys
import math
def extraction_atome(fichier_pdb, fichier_text):
    """ Lit un fichier pdb et renvoie les coordonées de tous les atomes
        impliqués dans la liaison peptidique (C,O,H,N).

    Parametres
    ----------
    fichier_pdb: fichier pdb contenant les coordonnées des atomes
    fichier_text: fichier de sortie qui contient que les coordonnées
                des atomes impliqués dans la liaison hydrogene

    Returns
    -------
    position: liste du numero des differents residus
    fichier_text

    """
    nom = ["C","O","N","H"]
    with open(fichier_pdb, "r") as pdb, open(fichier_text, "w") as text:
        for ligne in pdb:
            if ligne.startswith("ATOM") and ligne[12:16].strip() in nom:
                posi=ligne[22:27].strip()
                nam=ligne[12:16].strip()
                resid=ligne[17:20].strip()
                x = ligne[32:38]
                y = ligne[40:46]
                z = ligne[48:54]
                text.write("{}	{}	{}	{}	{}	{} \n".format(posi,resid,nam,x, y, z))
    return fichier_text


def dico_coord(file_text):
	dico = {}
	with open(file_text,"r") as t:
		for ligne in t:
			l = ligne.split("\t")
			posi = l[0]
			res = l[1]
			at = l[2]
			x, y, z = l[3], l[4], l[5].strip()
			if posi in dico:
				dico[posi].append((res, at, x,y,z))
				dico[posi] = [(res, at, x, y, z)]
	return dico

def distance(res1, res2):
    x1, y1, z1=float(res1[0]), float(res1[1]), float(res1[2])
    x2, y2, z2=float(res2[0]), float(res2[1]), float(res2[2])
    norme = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    distance=math.sqrt(norme)
    return distance


def calcul_energie(dico):
    cle=list(dico.keys())
    val=list(dico.values())
    q1 = -0.42
    q2 = 0.28
    f = 332
    liaison = {}
    for i in range(len(cle)-1):
        res1=dico[cle[i]]
        print(res1)
        C=res1[1][2:]
        O=res1[2][2:]
        #print(C,O)
        #xc,yc,zc=float(C[2]), float(C[3]), float(C[4])
        #xo,yo,zo=float(O[2]), float(O[3]), float(O[4])
        for j in range(i+1,len(cle)):
            res2 = dico[cle[j]]
            N = res2[0][2:]
            if len(res2)!=5:
            	H=(0,0,0)
            else:
            	H = res2[3][2:]
            ron=distance(O,N)
            rch=distance(C,H)
            roh=distance(O,H)
            rcn=distance(C,N)
            energie=q1*q2*f*(1/ron + 1/rch + 1/roh + 1/rcn)
            if energie < -0.5:
                liaison[(i,j)] = energie
    return liaison





if __name__ == '__main__':
	file_pdb = sys.argv[1]
	file_text = sys.argv[2]
	extraction_atome(file_pdb,file_text)
	dico=dico_coord(file_text)
	lh=calcul_energie(dico)
	print(lh)


























