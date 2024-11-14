.\bin\ymodel --scene tests\01_terrain\terrain.json --output outs\01_terrain\terrain.json --terrain object
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object
.\bin\ymodel --scene tests\03_hair1\hair1.json --output outs\03_hair1\hair1.json --hairbase object --hair hair
.\bin\ymodel --scene tests\03_hair2\hair2.json --output outs\03_hair2\hair2.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0
.\bin\ymodel --scene tests\03_hair3\hair3.json --output outs\03_hair3\hair3.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0.01
.\bin\ymodel --scene tests\03_hair4\hair4.json --output outs\03_hair4\hair4.json --hairbase object --hair hair --hairlen 0.02 --hairstr 0.001 --hairgrav 0.0005 --hairstep 8
.\bin\ymodel --scene tests\04_grass\grass.json --output outs\04_grass\grass.json --grassbase object --grass grass

.\bin\yscene render outs\01_terrain\terrain.json --output out\01_terrain.jpg --samples 256 --resolution  720
.\bin\yscene render outs\02_displacement\displacement.json --output out\02_displacement.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair1\hair1.json --output out\03_hair1.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair2\hair2.json --output out\03_hair2.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair3\hair3.json --output out\03_hair3.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair4\hair4.json --output out\03_hair4.jpg --samples 256 --resolution  720
.\bin\yscene render outs\04_grass\grass.json --output out\04_grass.jpg --samples 256 --resolution  720 --bounces 128

::precedural noises
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\05_voronoi\voronoi.json --displacement object --moded 2
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\06_voronoise\voronoise.json --displacement object --moded 3
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\07_cell\cell.json --displacement object --moded 1
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\08_smoothvoronoi\smoothvoronoi.json --displacement object --moded 4
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\09_realcell\realcell.json --displacement object --moded 5

.\bin\yscene render outs\05_voronoi\voronoi.json --output out\05_voronoi.jpg --samples 256 --resolution  720
.\bin\yscene render outs\06_voronoise\voronoise.json --output out\06_voronoise.jpg --samples 256 --resolution  720
.\bin\yscene render outs\07_cell\cell.json --output out\07_cell.jpg --samples 256 --resolution  720
.\bin\yscene render outs\08_smoothvoronoi\smoothvoronoi.json --output out\08_smoothvoronoi.jpg --samples 256 --resolution  720
.\bin\yscene render outs\09_realcell\realcell.json --output out\09_realcell.jpg --samples 256 --resolution  720


::Sample Elimination
.\bin\ymodel --scene tests\05_poisson\poisson.json --output outs\11_poisson_rabbit\poisson.json --hairbase object --hair hair --elimination --hairnum 1000 --hairlen 0.001 --hairstr 0
.\bin\yscene render outs\11_poisson_rabbit\poisson.json --output out\11_poisson_rabbit_many.jpg --samples 256 --resolution  720



::tree 
.\bin\ymodel --scene tests\06_trees\trees.json --output outs\12_trees\trees.json --trees object --steps 20
.\bin\yscene render outs\12_trees\trees.json --output out\12_trees_full_20_HD.jpg --samples 256 --resolution  720
.\bin\ymodel --scene tests\06_trees\trees.json --output outs\12_trees\trees.json --trees object --steps 100
.\bin\yscene render outs\12_trees\trees.json --output out\12_trees_full_100_HD.jpg --samples 256 --resolution  720
.\bin\ymodel --scene tests\06_trees\trees.json --output outs\12_trees\trees.json --trees object --steps 500
.\bin\yscene render outs\12_trees\trees.json --output out\12_trees_full_500_HD.jpg --samples 256 --resolution  720
.\bin\ymodel --scene tests\06_trees\trees.json --output outs\12_trees\trees.json --trees object --steps 1500
.\bin\yscene render outs\12_trees\trees.json --output out\12_trees_full_1500_HD.jpg --samples 256 --resolution  720
.\bin\ymodel --scene tests\06_trees\trees.json --output outs\12_trees\trees.json --trees object --steps 3500
.\bin\yscene render outs\12_trees\trees.json --output out\12_trees_full_2500_HD.jpg --samples 256 --resolution  720
.\bin\ymodel --scene tests\06_trees\trees.json --output outs\12_trees\trees.json --trees object --steps 10000
.\bin\yscene render outs\12_trees\trees.json --output out\12_trees_full_3500_HD.jpg --samples 256 --resolution  720

