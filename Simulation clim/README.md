![](schema%20projet%20clim.png)

->Main: Contient les fonctions et le programme se lance tout seul
   -Pour lancer le programme: py -i ./<chemin>/main.py (compilation standard en python)
   -Si vous souhaitez modifier des paramètres, entre les lignes 10 et 25 du main vous pouvez les modifiers manuellement.

 
->Le dossier FreeFem possède un fichier.edp pour crée les différents "cercle.msh" pour crée la canalisation.

->Le dossier msh possède le mesh du cercle:
	-"circle.msh" : nombre de triangle suffisant

->Le dossier "save" enregistre les résultats du fichier main. (Il est préférable de supprimer les images à l'interieur avant de relancé le programme si vous modifier le nombre d'itération).

![](save/sectionfinal.png)
