#!/usr/bin/env python3

#Solicitando o nome do arquivo ao usuário
arquivo = input("Digite o nome do arquivo FASTA: ")

#Criando os dicionários (vazios por agora)
fastaDict = {}  #Armazena os identificadores e suas sequências
ntCounts = {}   #Armazena, para cada identificador, sua contagem de nucleotídeos

#Abrindo o arquivo FASTA em modo leitura de texto
with open(arquivo, "rt") as fasta:
	for line in fasta:
		line = line.rstrip()  #Remove quebras de linha e espaços ao fim
		if line.startswith(">"):  #Identifica o cabeçalho de cada sequência
			identificador = line.split(" ")[0].replace(">","")  #Extrai o identificador da sequência, dividindo a linha pelo primeiro espaço. Também remove o caractere ">", para uma saída mais limpa.
			if identificador not in fastaDict.keys():
				#Caso o identificador ainda não seja uma chave no dicionário, inicia o valor com uma string vazia para acumular a sequência
				fastaDict[identificador] = ""
		else:
			#Adiciona o fragmento de sequência à sequência total correspondente ao identificador
			fastaDict[identificador] += line.upper()

#Percorrendo cada sequência armazenada no dicionário
for id in fastaDict.keys():
	nt_uniq = set(fastaDict[id])  #Obtém os nucleotídeos únicos presentes
	if id not in ntCounts.keys():
		ntCounts[id] = {}  #Cria um subdicionário para armazenar as contagens
	for nt in nt_uniq:
		if nt not in ntCounts[id]:
			ntCounts[id][nt] = fastaDict[id].count(nt)  #Conta quantas vezes cada nucleotídeo aparece na sequência

#Imprimindo o cabeçalho da tabela de saída
print("seqID\t\tA\tC\tG\tT")

#Percorrendo o dicionário de contagens e imprimindo os resultados
for id in ntCounts:
	#Usando 0 como valor padrão caso algum nucleotídeo esteja ausente (ex: sequência sem "G"), para evitar KeyError
	A = ntCounts[id]["A"] if "A" in ntCounts[id] else 0
	C = ntCounts[id]["C"] if "C" in ntCounts[id] else 0
	G = ntCounts[id]["G"] if "G" in ntCounts[id] else 0
	T = ntCounts[id]["T"] if "T" in ntCounts[id] else 0
	print(f"{id}\t{A}\t{C}\t{G}\t{T}")
