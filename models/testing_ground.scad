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