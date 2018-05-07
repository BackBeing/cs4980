`timescale 1ns / 1ps


module SoundLoopTest();

logic Areset, Clock, In;
logic [16:0] i;
logic [4:0] j;

// clock generation
parameter PERIOD = 10;
parameter pp = 17'b01111111111111111;

always
    begin Clock=1; #(PERIOD/2); Clock=0; #(PERIOD/2);
end


initial begin
    In = 0;
    Areset = 0; // hold reset signal for 5 cycles
    #(PERIOD*5);
    Areset = 1;
    
    // does In=0 disable the i --?
    #(PERIOD);
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
    #(PERIOD*5);
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
    #(PERIOD);
    
    // does the i-- works?
    In = 1;
    #(PERIOD); // s010
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
	#(PERIOD); // s100
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
	#(PERIOD); // s110
    assert (i == pp - 1) else $error("Expected: %d, Actual: %d", pp - 1 , i);
	#(PERIOD); // s010
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
	#(PERIOD); // s100
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
	#(PERIOD); // s110
    assert (i == pp - 1) else $error("Expected: %d, Actual: %d", pp - 1 , i);
    
    // does the Counter reset again?
    Areset = 0;
    #(PERIOD);
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
    
    // does the i-- again?
    Areset = 1;
    #(PERIOD); // s010
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
	#(PERIOD); // s100
    assert (i == pp) else $error("Expected: %d, Actual: %d", pp, i);
	#(PERIOD); // s110
    assert (i == pp - 1) else $error("Expected: %d, Actual: %d", pp - 1 , i);
end 

endmodule