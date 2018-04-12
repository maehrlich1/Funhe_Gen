"""
NGS Prep Clean-up
@author Opentrons
@date October 17th, 2017
@robot OT Pro, OT S
"""

from opentrons import containers, instruments, robot
import math

#These are max. speeds per motor! Do not increase further!
robot.head_speed(20000, z=6000, a=1200, b=1200)

"""
 Column A
"""
tip200_rack = containers.load('tiprack-200ul', 'A1')

H20 = containers.load('point', 'A2', 'H2O')

"""
 Column B
"""
trash = containers.load('point', 'B2', 'trash')

tip300_rack = containers.load('tiprack-200ul', 'B3')

"""
 Column C
"""
elution_plate = containers.load('96-PCR-flat', 'C1', 'elution plate')

EtOH = containers.load('point', 'C2', 'EtOH')

"""
 Column D
"""
gel_plate = containers.load('96-PCR-flat', 'D1', 'gel plate')

# trough with solutions
beads = containers.load('point', 'D2', 'beads')
# same trough with solutions
TE = containers.load('point', 'D2', 'TE')

mag_plate = containers.load('96-PCR-tall', 'D3', 'mag plate')

##############################################
# p200 (20 - 200 uL) (single)
p200single = instruments.Pipette(
    axis='b',
    name='p200single',
    max_volume=200,
    min_volume=20,
    channels=1,
    tip_racks=[tip200_rack])

p300multi = instruments.Pipette(
    axis='a',
    name='p300multi',
    max_volume=300,
    min_volume=50,
    channels=8,
    tip_racks=[tip300_rack])

mag_deck = instruments.Magbead(name='mag_deck')

# Helper Functions

def tip_wash(pipette, wash_vol, H2O_steps, EtOH_steps, lift_height):
	pipette.move_to((H20,(0,0,lift_height)))
	for _ in range(H2O_steps):
		pipette.move_to((H20,(0,0,lift_height)), strategy='direct')
		pipette.move_to(H20, strategy='direct')
		pipette.aspirate(wash_vol)
		pipette.move_to((H20,(0,0,lift_height)), strategy='direct')
		pipette.move_to((trash,(0,0,lift_height)), strategy='direct')
		pipette.move_to(trash, strategy='direct')
		pipette.dispense()
		pipette.move_to((trash,(0,0,lift_height)), strategy='direct')
	for _ in range(EtOH_steps):
		pipette.move_to((trash,(0,0,lift_height)), strategy='direct')
		pipette.move_to((EtOH,(0,0,lift_height)), strategy='direct')
		pipette.move_to(EtOH, strategy='direct')
		pipette.aspirate(wash_vol)
		pipette.move_to((EtOH,(0,0,lift_height)), strategy='direct')
		pipette.move_to((trash,(0,0,lift_height)), strategy='direct')
		pipette.move_to(trash, strategy='direct')
		pipette.dispense()

def EtOH_wash(pipette, plate, tiprack, vol, num_rows, discard_tip):
    p300multi.start_at_tip(tip300_rack.rows('1'))

    for row in range(num_rows):
        pipette.pick_up_tip()
        pipette.aspirate(vol, EtOH)
        pipette.dispense(vol, plate.rows(row))
        pipette.mix(5, vol)

        tip_wash(pipette, 200, 3, 2, 50)

        pipette.return_tip()

    mag_deck.engage()#.delay(minutes=5)

    p300multi.start_at_tip(tip300_rack.rows('1'))
    for row in range(num_rows):
        pipette.pick_up_tip()

        pipette.aspirate(vol, plate.rows(row))
        
        if discard_tip == False:
        	pipette.dispense(trash)

        	tip_wash(pipette, 200, 3, 2, 50)

        	pipette.return_tip()
        	
        else:
        	pipette.drop_tip(trash)
        	
#####################################################################        	
        	
def run_custom_protocol(TE_volume, number_of_samples):
    # all commands go in this function

    if(number_of_samples > 96):
        raise ValueError(
            "You can only input samples less than or equal to 96.")

    num_rows = math.ceil(number_of_samples/8)

    """
       Step 1 Add Magnetic Beads:
       -Get 8 tips from rack A
       -Go to bead trough and pipette up and do 5x with 100 ul to mix
       -Suck up 50 ul of beads
       -Dispense into PCR plate (magnet down),
       pipette up and down 5x with 100 ul to mix

       Step 2: Tip Wash
       -(Pick up 200 ul of water, expel into waste) x 3 times
       -(Pick up 200 ul of EtOH, expel into waste) x 2
       -Replace tips into box

    """
    
    # Step 1 and Step 2
    for row in range(num_rows):
        p300multi.pick_up_tip()
        p300multi.mix(5, 100, beads)
        p300multi.aspirate(50, beads)
        p300multi.dispense(50, mag_plate.rows(row))
        p300multi.mix(5, 100)

        tip_wash(p300multi, 200, 3, 2, 50)

        p300multi.return_tip()

    """
    Step 3: Bind and Beadwash
            -Bind DNA for 10 minutes
            -Raise magnent and wait 10 minutes
            -Remove 125 ul of elute from all wells
    Step 4: Wash tips
    """

    robot.home()
    #p300multi.delay(minutes=10)

    mag_deck.engage()#.delay(minutes=10)

    p300multi.start_at_tip(tip300_rack.rows('1'))
    for row in range(num_rows):
        p300multi.pick_up_tip()

        p300multi.aspirate(125, mag_plate.rows(row))
        p300multi.dispense(trash)

        tip_wash(p300multi, 200, 3, 2, 50)

        p300multi.return_tip()

    mag_deck.disengage()

    """
    Step 5: Wash all wells with EtOH
            -Add 200 EtOH to each well and pipette up and down
            5x
            -Wash tips
    Step 6: Lift magnent and wait five minutes
    Step 7: Remove 200 to waste
    """
    EtOH_wash(p300multi, mag_plate, tip300_rack, 200, num_rows, False)

    """
    Step 8: Disengage mag deck
    Step 9: Wash all wells with 200 EtOH
    Step 10: Disengage mag deck
    Step 11: Wash all wells with 200 EtOH
    Step 12: Disengage mag deck
    """
    mag_deck.disengage()

    EtOH_wash(p300multi, mag_plate, tip300_rack, 200, num_rows, False)

    mag_deck.disengage()

    EtOH_wash(p300multi, mag_plate, tip300_rack, 200, num_rows, True)

    mag_deck.disengage()

    robot.home()
    #p300multi.delay(minutes=5)

    """
    Step 13: Elute DNA from all rows
             -Resuspend in TE
             -Sit 5 minutes
             -Magnetize
             -Put elute in clean plate and gel plate
    """
    p200single.start_at_tip(tip200_rack.wells('A1'))
    for row in range(num_rows):
        for well in mag_plate.rows(row):
            p200single.pick_up_tip()

            p200single.aspirate(TE_volume+5, TE)
            p200single.dispense(TE_volume+5, well)
            p200single.mix(5, TE_volume+5)

            tip_wash(p200single, 200, 3, 2, 50)

            p200single.return_tip()

    #p200single.delay(minutes=5)

    mag_deck.engage()#.delay(minutes=5)

    p200single.start_at_tip(tip200_rack.wells('A1'))

    for i in range(num_rows):
        for well1, well2, well3 in zip(mag_plate.rows(i),
                                       elution_plate.rows(i),
                                       gel_plate.rows(i)):
            p200single.pick_up_tip()

            p200single.aspirate(TE_volume+5, well1)
            p200single.dispense(TE_volume, well2)

            p200single.dispense(5, well3)

            p200single.drop_tip(trash)

    mag_deck.disengage()
    
    robot.home()

run_custom_protocol(**{'TE_volume': 20.0, 'number_of_samples': 8})