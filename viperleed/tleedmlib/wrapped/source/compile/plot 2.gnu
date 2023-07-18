set terminal pdfcairo enhanced font "Arial,12"
set encoding iso_8859_1
set xrange [ 71:799]
set yrange [0:1499]
set size 1.000,1.000
set style data lines
set xtics  50
set noytics
set xlabel "Energy [eV]"
set noxzeroaxis

set output "pic_001.pdf"
set nolabel
set label "*( 1  0)*           \nR = 0.3693" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($  1>0?$  1:1/0) :($  1>0?$  2:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($  1>0?$  1:1/0) :($  1>0?$  2:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_002.pdf"
set nolabel
set label "*(-1 -1)*           \nR = 0.2038" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($  3>0?$  3:1/0) :($  3>0?$  4:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($  3>0?$  3:1/0) :($  3>0?$  4:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_003.pdf"
set nolabel
set label "*( 1  1)*           \nR = 0.1002" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($  5>0?$  5:1/0) :($  5>0?$  6:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($  5>0?$  5:1/0) :($  5>0?$  6:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_004.pdf"
set nolabel
set label "*( 0 -2)*           \nR = 0.1051" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($  7>0?$  7:1/0) :($  7>0?$  8:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($  7>0?$  7:1/0) :($  7>0?$  8:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_005.pdf"
set nolabel
set label "*( 0  2)*           \nR = 0.2062" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($  9>0?$  9:1/0) :($  9>0?$ 10:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($  9>0?$  9:1/0) :($  9>0?$ 10:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_006.pdf"
set nolabel
set label "*( 2  0)*           \nR = 0.1407" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 11>0?$ 11:1/0) :($ 11>0?$ 12:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 11>0?$ 11:1/0) :($ 11>0?$ 12:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_007.pdf"
set nolabel
set label "*(-1 -2)*           \nR = 0.1956" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 13>0?$ 13:1/0) :($ 13>0?$ 14:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 13>0?$ 13:1/0) :($ 13>0?$ 14:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_008.pdf"
set nolabel
set label "*( 1  2)*           \nR = 0.3957" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 15>0?$ 15:1/0) :($ 15>0?$ 16:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 15>0?$ 15:1/0) :($ 15>0?$ 16:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_009.pdf"
set nolabel
set label "*( 2  1)*           \nR = 0.1248" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 17>0?$ 17:1/0) :($ 17>0?$ 18:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 17>0?$ 17:1/0) :($ 17>0?$ 18:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_010.pdf"
set nolabel
set label "*( 2 -1)*           \nR = 0.2067" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 19>0?$ 19:1/0) :($ 19>0?$ 20:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 19>0?$ 19:1/0) :($ 19>0?$ 20:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_011.pdf"
set nolabel
set label "*( 2 -2)*           \nR = 0.1481" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 21>0?$ 21:1/0) :($ 21>0?$ 22:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 21>0?$ 21:1/0) :($ 21>0?$ 22:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_012.pdf"
set nolabel
set label "*( 2  2)*           \nR = 0.2597" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 23>0?$ 23:1/0) :($ 23>0?$ 24:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 23>0?$ 23:1/0) :($ 23>0?$ 24:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_013.pdf"
set nolabel
set label "*( 1 -3)*           \nR = 0.1845" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 25>0?$ 25:1/0) :($ 25>0?$ 26:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 25>0?$ 25:1/0) :($ 25>0?$ 26:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_014.pdf"
set nolabel
set label "*( 1  3)*           \nR = 0.2749" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 27>0?$ 27:1/0) :($ 27>0?$ 28:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 27>0?$ 27:1/0) :($ 27>0?$ 28:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_015.pdf"
set nolabel
set label "*( 3  0)*           \nR = 0.2362" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 29>0?$ 29:1/0) :($ 29>0?$ 30:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 29>0?$ 29:1/0) :($ 29>0?$ 30:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_016.pdf"
set nolabel
set label "*( 3  1)*           \nR = 0.2661" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 31>0?$ 31:1/0) :($ 31>0?$ 32:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 31>0?$ 31:1/0) :($ 31>0?$ 32:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_017.pdf"
set nolabel
set label "*( 3 -1)*           \nR = 0.3155" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 33>0?$ 33:1/0) :($ 33>0?$ 34:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 33>0?$ 33:1/0) :($ 33>0?$ 34:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_018.pdf"
set nolabel
set label "*(-2 -3)*           \nR = 0.3363" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 35>0?$ 35:1/0) :($ 35>0?$ 36:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 35>0?$ 35:1/0) :($ 35>0?$ 36:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_019.pdf"
set nolabel
set label "*( 2  3)*           \nR = 0.4154" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 37>0?$ 37:1/0) :($ 37>0?$ 38:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 37>0?$ 37:1/0) :($ 37>0?$ 38:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_020.pdf"
set nolabel
set label "*( 3  2)*           \nR = 0.2750" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 39>0?$ 39:1/0) :($ 39>0?$ 40:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 39>0?$ 39:1/0) :($ 39>0?$ 40:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_021.pdf"
set nolabel
set label "*( 3 -2)*           \nR = 0.3806" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 41>0?$ 41:1/0) :($ 41>0?$ 42:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 41>0?$ 41:1/0) :($ 41>0?$ 42:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_022.pdf"
set nolabel
set label "*( 0 -4)*           \nR = 0.3709" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 43>0?$ 43:1/0) :($ 43>0?$ 44:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 43>0?$ 43:1/0) :($ 43>0?$ 44:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_023.pdf"
set nolabel
set label "*( 0  4)*           \nR = 0.4518" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 45>0?$ 45:1/0) :($ 45>0?$ 46:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 45>0?$ 45:1/0) :($ 45>0?$ 46:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_024.pdf"
set nolabel
set label "*( 1 -4)*           \nR = 0.2681" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 47>0?$ 47:1/0) :($ 47>0?$ 48:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 47>0?$ 47:1/0) :($ 47>0?$ 48:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_025.pdf"
set nolabel
set label "*( 1  4)*           \nR = 0.6442" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 49>0?$ 49:1/0) :($ 49>0?$ 50:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 49>0?$ 49:1/0) :($ 49>0?$ 50:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_026.pdf"
set nolabel
set label "*( 4  0)*           \nR = 0.2327" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 51>0?$ 51:1/0) :($ 51>0?$ 52:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 51>0?$ 51:1/0) :($ 51>0?$ 52:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_027.pdf"
set nolabel
set label "*( 3  3)*           \nR = 0.2839" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 53>0?$ 53:1/0) :($ 53>0?$ 54:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 53>0?$ 53:1/0) :($ 53>0?$ 54:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_028.pdf"
set nolabel
set label "*( 3 -3)*           \nR = 0.3071" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 55>0?$ 55:1/0) :($ 55>0?$ 56:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 55>0?$ 55:1/0) :($ 55>0?$ 56:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_029.pdf"
set nolabel
set label "*( 4 -1)*           \nR = 0.4120" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 57>0?$ 57:1/0) :($ 57>0?$ 58:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 57>0?$ 57:1/0) :($ 57>0?$ 58:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_030.pdf"
set nolabel
set label "*( 4  1)*           \nR = 0.4857" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 59>0?$ 59:1/0) :($ 59>0?$ 60:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 59>0?$ 59:1/0) :($ 59>0?$ 60:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_031.pdf"
set nolabel
set label "*(-2 -4)*           \nR = 0.2522" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 61>0?$ 61:1/0) :($ 61>0?$ 62:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 61>0?$ 61:1/0) :($ 61>0?$ 62:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_032.pdf"
set nolabel
set label "*( 2  4)*           \nR = 0.4216" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 63>0?$ 63:1/0) :($ 63>0?$ 64:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 63>0?$ 63:1/0) :($ 63>0?$ 64:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_033.pdf"
set nolabel
set label "*( 4  2)*           \nR = 0.4203" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 65>0?$ 65:1/0) :($ 65>0?$ 66:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 65>0?$ 65:1/0) :($ 65>0?$ 66:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_034.pdf"
set nolabel
set label "*( 4 -2)*           \nR = 0.3531" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 67>0?$ 67:1/0) :($ 67>0?$ 68:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 67>0?$ 67:1/0) :($ 67>0?$ 68:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_035.pdf"
set nolabel
set label "*( 1 -5)*           \nR = 0.3997" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 69>0?$ 69:1/0) :($ 69>0?$ 70:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 69>0?$ 69:1/0) :($ 69>0?$ 70:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_036.pdf"
set nolabel
set label "*( 1  5)*           \nR = 0.3275" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 71>0?$ 71:1/0) :($ 71>0?$ 72:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 71>0?$ 71:1/0) :($ 71>0?$ 72:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_037.pdf"
set nolabel
set label "*(-3 -4)*           \nR = 0.8673" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 73>0?$ 73:1/0) :($ 73>0?$ 74:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 73>0?$ 73:1/0) :($ 73>0?$ 74:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_038.pdf"
set nolabel
set label "*( 3  4)*           \nR = 0.4087" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 75>0?$ 75:1/0) :($ 75>0?$ 76:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 75>0?$ 75:1/0) :($ 75>0?$ 76:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_039.pdf"
set nolabel
set label "*( 4  3)*           \nR = 0.4449" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 77>0?$ 77:1/0) :($ 77>0?$ 78:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 77>0?$ 77:1/0) :($ 77>0?$ 78:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_040.pdf"
set nolabel
set label "*( 4 -3)*           \nR = 0.5570" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 79>0?$ 79:1/0) :($ 79>0?$ 80:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 79>0?$ 79:1/0) :($ 79>0?$ 80:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_041.pdf"
set nolabel
set label "*( 5  0)*           \nR = 0.4042" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 81>0?$ 81:1/0) :($ 81>0?$ 82:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 81>0?$ 81:1/0) :($ 81>0?$ 82:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_042.pdf"
set nolabel
set label "*( 2 -5)*           \nR = 0.4109" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 83>0?$ 83:1/0) :($ 83>0?$ 84:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 83>0?$ 83:1/0) :($ 83>0?$ 84:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_043.pdf"
set nolabel
set label "*( 2  5)*           \nR = 0.4702" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 85>0?$ 85:1/0) :($ 85>0?$ 86:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 85>0?$ 85:1/0) :($ 85>0?$ 86:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_044.pdf"
set nolabel
set label "*( 5  1)*           \nR = 0.4945" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 87>0?$ 87:1/0) :($ 87>0?$ 88:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 87>0?$ 87:1/0) :($ 87>0?$ 88:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_045.pdf"
set nolabel
set label "*( 5 -1)*           \nR = 0.5288" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 89>0?$ 89:1/0) :($ 89>0?$ 90:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 89>0?$ 89:1/0) :($ 89>0?$ 90:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_046.pdf"
set nolabel
set label "*( 5 -2)*           \nR = 0.5541" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 91>0?$ 91:1/0) :($ 91>0?$ 92:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 91>0?$ 91:1/0) :($ 91>0?$ 92:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_047.pdf"
set nolabel
set label "*( 5  2)*           \nR = 0.6924" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 93>0?$ 93:1/0) :($ 93>0?$ 94:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 93>0?$ 93:1/0) :($ 93>0?$ 94:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_048.pdf"
set nolabel
set label "*( 4 -4)*           \nR = 0.8519" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 95>0?$ 95:1/0) :($ 95>0?$ 96:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 95>0?$ 95:1/0) :($ 95>0?$ 96:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_049.pdf"
set nolabel
set label "*( 4  4)*           \nR = 0.2266" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 97>0?$ 97:1/0) :($ 97>0?$ 98:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 97>0?$ 97:1/0) :($ 97>0?$ 98:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_050.pdf"
set nolabel
set label "*(-3 -5)*           \nR = 0.5075" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($ 99>0?$ 99:1/0) :($ 99>0?$100:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($ 99>0?$ 99:1/0) :($ 99>0?$100:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_051.pdf"
set nolabel
set label "*( 3  5)*           \nR = 0.4724" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($101>0?$101:1/0) :($101>0?$102:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($101>0?$101:1/0) :($101>0?$102:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_052.pdf"
set nolabel
set label "*( 0 -6)*           \nR = 0.4576" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($103>0?$103:1/0) :($103>0?$104:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($103>0?$103:1/0) :($103>0?$104:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_053.pdf"
set nolabel
set label "*( 0  6)*           \nR = 0.5167" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($105>0?$105:1/0) :($105>0?$106:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($105>0?$105:1/0) :($105>0?$106:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_054.pdf"
set nolabel
set label "*( 5 -3)*           \nR = 0.2246" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($107>0?$107:1/0) :($107>0?$108:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($107>0?$107:1/0) :($107>0?$108:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_055.pdf"
set nolabel
set label "*( 5  3)*           \nR = 0.5595" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($109>0?$109:1/0) :($109>0?$110:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($109>0?$109:1/0) :($109>0?$110:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_056.pdf"
set nolabel
set label "*( 2 -6)*           \nR = 0.5643" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($111>0?$111:1/0) :($111>0?$112:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($111>0?$111:1/0) :($111>0?$112:1/0) title "(theory)" lt rgb "blue" lw 2

set output "pic_057.pdf"
set nolabel
set label "*( 2  6)*           \nR = 0.6365" at 435,1199 center font "Helvetica,12"
plot 'exp.column' using ($113>0?$113:1/0) :($113>0?$114:1/0) title "(exp.)" lt rgb "red" lw 2,\
 'theo.column' using ($113>0?$113:1/0) :($113>0?$114:1/0) title "(theory)" lt rgb "blue" lw 2

