/*intersection() {    
    difference() {
        translate([0, 0, 0])
            sphere(r=8.5);
        translate([0, 5, 0])
            sphere(r=8.5);
        translate([0, -5, 0])
            sphere(r=8.5);
    }
}*/

/*intersection() {    
    difference() {
        translate([0, 0, 0])
            sphere(r=8.5);
        translate([0, 5, 5])
            sphere(r=8.5);
        translate([0, -5, 0])
            sphere(r=8.5);
    }
}*/

/*
difference() {
    
rotate([0, 20, 0]) {
    difference() {
        translate([0, 0, 0])
            sphere(r=8.5);
        translate([0, 0, 8])
            sphere(r=8.5);
        translate([0, 0, -8])
            sphere(r=8.5);
       translate([15, 0, 3])
            sphere(r=8.5);
        translate([0, 15, 3])
            sphere(r=8.5);
        
    }

}
    translate([-15, -15, -1.4])
            cube([30, 30, 20]);
}
*/


intersection() {
union() {   
translate([-5, -5, 0])
    cube([10, 10, 10]);
translate([3, -5, 9])
    cube([10, 10, 10]);
translate([-13, -5, 6])
    cube([10, 10, 10]);
translate([-5, -5, 15])
    cube([10, 10, 10]);
}

//translate([-50, -50, 0])
//cube([100, 100, 1]);

//translate([-50, -50, 5])
//cube([100, 100, 2]);

//translate([-50, -50, 8])
//cube([100, 100, 2]);

//translate([-50, -50, 9])
//cube([100, 100, 2]);

//translate([-50, -50, 14])
//cube([100, 100, 2]);

//translate([-50, -50, 14])
//cube([100, 100, 2]);

//translate([-50, -50, 15])
//cube([100, 100, 2]);

//translate([-50, -50, 18])
//cube([100, 100, 2]);

//translate([-50, -50, 24])
//cube([100, 100, 2]);
}
