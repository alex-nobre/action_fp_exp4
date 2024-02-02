# -*- coding: utf-8 -*-
"""
Checks for possible delays between the action trial start keypress and the presentation of the warning signal.
"""

# Import libraries
import pandas as pd
import glob
import os


# Instruction texts for comparisons
action_start = "instrucoes_practica_action: text = 'Nesta primeira parte da tarefa, a cruz preta aparecerá na tela e permanecerá lá até você pressionar a tecla 1 na caixa de respostas para fazê-la mudar de cor. Utilize sua mão não-dominante para pressionar a tecla 1 (p. ex., se você é destro, utilize a mão esquerda para pressionar a tecla 1). Só então a cruz mudará de preta para verde. Nesse momento, você deverá se preparar para a apresentação do estímulo.\\n\\nDepois de responder ao estímulo, aguarde até que a cruz preta apareça novamente antes de apertar a tecla 1 para fazê-la mudar de cor novamente.\\n\\nAntes de iniciarmos, vamos praticar um pouco. Pressione a tecla 1 para iniciar a fase de prática.'"
external_start = "instrucoes_pratica_external: text = 'Nesta primeira parte da tarefa, a cruz mudará de cor automaticamente. Você deverá ficar atento ao momento em que a mudança de cor acontecer e se preparar para a apresentação do estímulo.\\n\\nAntes de iniciarmos, vamos praticar um pouco. Pressione a tecla 1 na caixa de respostas para iniciar a fase de prática.'"

action_middle = "instrucoes_practica_action: text = 'Nos próximos blocos, a tarefa será um pouco diferente. Ao invés de esperar a cruz mudar de cor sozinha, você deverá pressionar a tecla 1 na caixa de respostas para fazê-la mudar de cor (de preta para verde). Utilize sua mão não-dominante para pressionar a tecla 1 (p. ex., se você é destro, utilize a mão esquerda para pressionar a tecla 1).\\n\\nDepois de responder ao estímulo, aguarde até que a cruz preta apareça novamente antes apertar a tecla 1 para fazê-la mudar de cor novamente.\\n\\nNos outros aspectos, a tarefa será igual a antes.\\n\\nPrimeiro, pratique algumas vezes para ter certeza de que compreendeu a mudança. Pressione a tecla 1 na caixa de respostas para iniciar a fase de prática.'"
external_middle = "instrucoes_pratica_external: text = 'Nos próximos blocos, a tarefa será um pouco diferente. Ao invés de pressionar a tecla 1 para fazer a cruz mudar de cor, ela mudará de cor sozinha. Você deverá ficar atento ao momento em que a mudança de cor acontecer.\\n\\nNos outros aspectos, a tarefa será igual a antes.\\n\\nPrimeiro, pratique algumas vezes para ter certeza de que compreendeu a mudança. Pressione a tecla 1 para iniciar a fase de prática.'"

#practice_start = 'Antes de iniciarmos, vamos praticar um pouco. Pressione a tecla 1 na caixa de respostas para iniciar a fase de prática.'
#practice_middle = 'Primeiro, pratique algumas vezes para ter certeza de que compreendeu a mudança. Pressione a tecla 1 na caixa de respostas para iniciar a fase de prática.'

practice_end_once = 'A fase de prática terminou. Agora iremos iniciar a tarefa. Ela será igual à prática, porém mais longa.'
practice_end_repeat = 'Vamos repetir a fase de prática para ter certeza de que você compreendeu as instruções.'

# Path to log files
filesPath = '.\Data\Log_files'

# Get list of log files
FileList=glob.glob(filesPath + '/*.log')
FileList.sort()

nFiles=int(len(FileList))

# Delay data file for all subs
delayDataAll = pd.DataFrame(columns = ['participant', 'condition',  'trial', 'delay'])

for iFile, fileName in enumerate(FileList):
    
    # read file in correct format 
    logData = open(fileName).read()
    
    # Split to read line by line
    lines = logData.split('\n')
    numberLines = len(lines)
    
    # Get participant ID from file name
    participant = fileName[17:20]
    
    # Set practice to False (changes program behavior when set to True)
    practice = False
    
    # Create placeholder condition variable
    condition = None
    
    # Trial counter
    trialCount = 1
    
    # Begin outside trial
    trialStart = False
    
    # Dataset to save delay info for subject
    subDelayData = pd.DataFrame(columns = ['participant', 'condition',  'trial', 'delay'])
    subFileName = participant + '_delayData' + '.csv'
    
    print(participant)
    
    # Read line by line
    for line in range(numberLines - 1):
      # Split line into columns (only the third columns is of interest)
      columns = lines[line].split('\t')
      
      # Search for test phase start on third column
      if practice == False:
        if len(columns) == 3:
          # Check for practice start and set practice to True if started; does nothing otherwise
          for thisText in [action_start, action_middle, external_start, external_middle]:
            if columns[2].find(thisText) >= 0:
              practice = True # does not look for responses
              if columns[2] in [external_start, external_middle]:
                condition = 'external'
              elif columns[2] in [action_start, action_middle]:
                condition = 'action'
            
      if practice == False: # check again before running action checks 
        # Run through trial lines
        if len(columns) == 3:
          # For action condition, look up start of keydown time for delay start time and fixation 2 draw start for delay end time
          if condition == "action":
            if trialStart == False:
              if columns[2] in ['Fixation_action: autoDraw = True']: # set trialstart here to avoid confusion with response keydown
                trialStart = True # allows look up of delay end time
            elif trialStart == True:
              if columns[2] in ['Keypress: 1']:
                delayStart = columns[0]
              elif columns[2] in ['Fixation_2: autoDraw = True']:
                delayEnd = columns[0]
                delay = float(delayEnd) - float(delayStart)
                trialStart = False # allows look up of delay start time
                
                # Add to dataset
                trialDelayData = [participant, condition, trialCount, delay]
                subDelayData.loc[len(subDelayData.index)] = trialDelayData
          
                # update trial counter
                trialCount += 1
                
          # For external condition, look up start of external fixation for delay start time and fixation 2 draw start for delay end time
          if condition == 'external':
            if trialStart == False:
              if columns[2] in ['Fixation_ext: autoDraw = False']:
                trialStart = True # allows look up of trial end time
                delayStart = columns[0]
            elif trialStart == True:
              if columns[2] in ['Fixation_2: autoDraw = True']:
                delayEnd = columns[0]
                delay = float(delayEnd) - float(delayStart)
                trialStart = False # allows look up of trial start time
          
                # Add to dataset
                trialDelayData = [participant, condition, trialCount, delay]
                subDelayData.loc[len(subDelayData.index)] = trialDelayData
          
                # update trial counter
                trialCount += 1
                
      # Check if practice has ended and determine condition
      if practice == True:
        if len(columns) == 3:
          # Check that practice was not repeated
          if columns[2].find(practice_end_once) > 0:
            practice = False
          elif columns[2].find(practice_end_repeat) > 0:
            practice = False

    subDelayData.to_csv('./Analysis/Delay_data/' + subFileName)
    
    # Add to global dataset
    delayDataAll = pd.concat([delayDataAll, subDelayData], axis = 0)
    

# Save global dataset
delayDataAll.to_csv('./Analysis/'+'delayDataAll.csv')
