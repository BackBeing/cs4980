`timescale 1ns / 1ps


module CounterTest();

logic Areset, Clock, In;
logic [16:0] i;
logic [4:0] j;

// clock generation
parameter PERIOD = 10;
parameter pp = 17'b01111111111111111;

always
    begin Clock=1; #(PERIOD/2); Clock=0; #(PERIOD/2);
end

// instantiate a 5-bit version of Counter
Counter #(4) dut(.Areset(Areset), .Clock(Clock), .In(In), .i(i), .j(j));

initial begin
    In = 0;
    Areset = 0; // hold reset signal for 5 cycles
    #(PERIOD*5);
    Areset = 1;
    
    // does In=0 disable the i --?
    #(PERIOD);
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
    #(PERIOD*5);
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
    #(PERIOD);
    
    // does the i-- works?
    In = 1;
    #(PERIOD); // s010
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
	#(PERIOD); // s100
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
	#(PERIOD); // s110
    assert (i == pp - 1) else $error("Expected: %d, Actual: %d", pp - 1 , Count);
	#(PERIOD); // s010
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
	#(PERIOD); // s100
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
	#(PERIOD); // s110
    assert (i == pp - 1) else $error("Expected: %d, Actual: %d", pp - 1 , Count);
    
    // does the Counter reset again?
    Areset = 0;
    #(PERIOD);
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
    
    // does the i-- again?
    Areset = 1;
    #(PERIOD); // s010
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
	#(PERIOD); // s100
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, Count);
	#(PERIOD); // s110
    assert (i == pp - 1) else $error("Expected: %d, Actual: %d", pp - 1 , Count);
end 

endmodule