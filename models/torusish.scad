intersection() {    
    difference() {
        rotate([90, 0, 0]) {
            translate([0, 0, -5]) {
                cylinder(10, 10, 10);
            }
        }
        translate([0, -5, 0])
            sphere(r=8.5);
        translate([0, 5, 0])
            sphere(r=8.5);
    }
    //sphere(r=10.0);
}