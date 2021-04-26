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

intersection() {    
    difference() {
        translate([0, 0, 0])
            sphere(r=8.5);
        translate([0, 0, 8])
            sphere(r=8.5);
        translate([0, 0, -8])
            sphere(r=8.5);
        translate([0, 3, 0])
            sphere(r=8.5);
    }
}