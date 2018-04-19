`timescale 1ns / 1ps


// This module is broken, you must debug it
module SoundLoop(
    input In,
    output Out,
    input Clock,
    input Areset
    );
    
    // state encodings
    parameter s00 = 2'b00;
    parameter s01 = 2'b01;
    parameter s10 = 2'b10;
    parameter s11 = 2'b11;
    
    // 2-bit register holding the current state of the FSM
    logic [1:0] state;
    
    // 2-bit combinational logic result
    logic [1:0] nextState;
    
    always_ff@(posedge Clock) begin
        if (Areset == 1'b0) begin
            state <= s00;
        end else begin
            state <= nextState;
        end
    end
    
    // Output logic
    assign Out = (state == s11) ? 1'b1 : 1'b0;
    
    // Next state logic
    always_comb begin
        case(state)
        s00 : if (In == 1'b0) nextState = s00;
            else nextState = s10;
        s01 : if (In == 1'b0) nextState = s10;
            else nextState = s11;
        s10 : if (In == 1'b0) nextState = s00;
            else nextState = s01;
        s11 : if (In == 1'b0) nextState = s01;
            else nextState = s11;
        endcase
    end
endmodule
