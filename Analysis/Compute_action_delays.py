# -*- coding: utf-8 -*-
"""
Checks for possible delays between the action trial start keypress and the presentation of the warning signal.
"""

# Import libraries
import pandas as pd
import glob
import os


# Instruction texts for comparisons
action_start = 'instrucoes_practica_action: text = Nesta primeira parte da tarefa, a cruz preta aparecerá na tela e permanecerá lá até você pressionar a tecla 1 na caixa de respostas para fazê-la mudar de cor. Utilize sua mão não-dominante para pressionar a tecla 1 (p. ex., se você é destro, utilize a mão esquerda para pressionar a tecla 1). Só então a cruz mudará de preta para verde. Nesse momento, você deverá se preparar para a apresentação do estímulo.'
external_start = 'instrucoes_pratica_external: text = Nesta primeira parte da tarefa, a cruz mudará de cor automaticamente. Você deverá ficar atento ao momento em que a mudança de cor acontecer e se preparar para a apresentação do estímulo.'

action_middle = 'instrucoes_practica_action: text = Nos próximos blocos, a tarefa será um pouco diferente. Ao invés de esperar a cruz mudar de cor sozinha, você deverá pressionar a tecla 1 na caixa de respostas para fazê-la mudar de cor (de preta para verde). Utilize sua mão não-dominante para pressionar a tecla 1 (p. ex., se você é destro, utilize a mão esquerda para pressionar a tecla 1).'
external_middle = 'instrucoes_pratica_external: text = Nos próximos blocos, a tarefa será um pouco diferente. Ao invés de pressionar a tecla 1 para fazer a cruz mudar de cor, ela mudará de cor sozinha. Você deverá ficar atento ao momento em que a mudança de cor acontecer.'

practice_start = 'Antes de iniciarmos, vamos praticar um pouco. Pressione a tecla 1 na caixa de respostas para iniciar a fase de prática.'
practice_middle = 'Primeiro, pratique algumas vezes para ter certeza de que compreendeu a mudança. Pressione a tecla 1 na caixa de respostas para iniciar a fase de prática.'

practice_end_once = 'A fase de prática terminou. Agora iremos iniciar a tarefa. Ela será igual à prática, porém mais longa.'
practice_end_repeat = 'Vamos repetir a fase de prática para ter certeza de que você compreendeu as instruções.'

# Path to log files
filesPath = '.\Data\Log_files'

# Get list of log files
FileList=glob.glob(filesPath + '/*.log')
FileList.sort()

nFiles=int(len(FileList))

# Delay data file for all subs
delayDataAll = pd.DataFrame(columns = ['ID', 'condition',  'trial', 'delay'])

for iFile, fileName in enumerate(FileList):
    
    # read file in correct format 
    logData = open(fileName).read()
    
    # Split to read line by line
    lines = logData.split('\n')
    numberLines = len(lines)
    
    # Get participant ID from file name
    ID = fileName[17:20]
    
    # Set practice to False (changes program behavior when set to True)
    practice = False
    
    # Create placeholder condition variable
    condition = None
    
    # Trial counter
    trialCount = 1
    
    # Begin outside trial
    trialStart = False
    
    # Dataset to save delay info for subject
    subDelayData = pd.DataFrame(columns = ['ID', 'condition',  'trial', 'delay'])
    subFileName = ID + '_delayData' + '.csv'
    
    print(ID)
    
    # Read line by line
    for line in range(numberLines - 1):
      # Split line into columns (only the third columns is of interest)
      columns = lines[line].split('\t')
      
      # Perform comparisons on third column
      if practice == False:
        if len(columns) == 1:
          # Check for practice start and set practice to True if started; does nothing otherwise
          if columns[0] in [practice_start, practice_middle]:
            practice = True # does not look for responses
            
        
        # Run through trial lines
        if len(columns) == 3:
          # For action condition, look up start of keydown time for delay start time and fixation 2 draw start for delay end time
          if condition == "action":
            if trialStart == False:
              if columns[2] in ['Fixation_action: autoDraw = true']: # set trialstart here to avoid confusion with response keydown
                trialStart = True # allows look up of delay end time
                #print('EXECUTED')
            elif trialStart == True:
              if columns[2] in ['Keypress: 1']:
                delayStart = columns[0]
              elif columns[2] in ['Fixation_2: autoDraw = true']:
                delayEnd = columns[0]
                delay = float(delayEnd) - float(delayStart)
                trialStart = False # allows look up of delay start time
                
                # Add to dataset
                trialDelayData = [ID, condition, trialCount, delay]
                subDelayData.loc[len(subDelayData.index)] = trialDelayData
          
                # update trial counter
                #print('external ' + str(trialCount))
                trialCount += 1
                
          # For external condition, look up start of external fixation for delay start time and fixation 2 draw start for delay end time
          if condition == 'external':
            if trialStart == False:
              if columns[2] in ['Fixation_ext: autoDraw = False']:
                trialStart = True # allows look up of trial end time
                delayStart = columns[0]
            elif trialStart == True:
              if columns[2] in ['Fixation_2: autoDraw = true']:
                delayEnd = columns[0]
                delay = float(delayEnd) - float(delayStart)
                trialStart = False # allows look up of trial start time
          
                # Add to dataset
                trialDelayData = [ID, condition, trialCount, delay]
                subDelayData.loc[len(subDelayData.index)] = trialDelayData
          
                # update trial counter
                #print('action ' + str(trialCount))
                trialCount += 1
                
        #if len(columns) == 3:
          # Determine condition based on text
          if columns[2] in [external_start, external_middle]:
            condition = 'external'
            print('yes')
          elif columns[2] in [action_start, action_middle]:
            condition = 'action'
            print('yes')
        
      # Check if practice has ended and determine condition
      if practice == True:
        #print('practiceTrue')
        if len(columns) == 3:
          # Check that practice was not repeated
          if columns[2] in [practice_end_once, practice_end_repeat]:
          #if columns[0] == 'A fase de prática terminou. Agora iremos iniciar a tarefa. Ela será igual à prática, porém mais longa.':
            print('match')
            practice = False
    

    subDelayData.to_csv('./Analysis/Delay_data/' + subFileName)
    
    # Add to global dataset
    delayDataAll = pd.concat([delayDataAll, subDelayData], axis = 0)
    

# Save global dataset
delayDataAll.to_csv('./Analysis/'+'delayDataAll.csv')
