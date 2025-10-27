#!/usr/bin/env python3

#Solicitando o nome do arquivo ao usuário
arquivo = input("Digite o nome do arquivo FASTA: ")

#Criando os dicionários (vazios por agora)
fastaDict = {}  #Armazena os identificadores e suas sequências

#Abrindo o arquivo FASTA (leitura) e o arquivo de saída (escrita)
with open(arquivo, "rt") as fasta, open("Python_08.codons-3frames.nt", "w") as arquivo_saida:
	for line in fasta:
		line = line.rstrip()  #Remove quebras de linha e espaços ao fim

		if line.startswith(">"):  #Identifica o cabeçalho de cada sequência
			identificador = line.split(" ")[0].replace(">","")  #Extrai o identificador da sequência, dividindo a linha pelo primeiro espaço. Também remove o caractere ">", para uma saída mais limpa.
			if identificador not in fastaDict.keys():
				#Caso o identificador ainda não seja uma chave no dicionário, inicia o valor com uma string vazia para acumular a sequência
				fastaDict[identificador] = ""

		else:
			#Adiciona o fragmento de sequência à sequência total correspondente ao identificador
			#Certifica que a sequência está em maiúsculas
			fastaDict[identificador] += line.upper()


	#Nesse ponto, o arquivo FASTA foi lido e o dicionário fastaDict está completo
	#Agora, iteramos pelo dicionário para processar cada sequência e escrever no arquivo de saída
	for identificador, sequencia in fastaDict.items():

		#Criando um novo loop para iterar pelos distintos quadros de leitura
		for quadro in range(1,4):
			inicio_index = quadro-1
			arquivo_saida.write(f"{identificador}-frame-{quadro}-codons\n")  #Escreve o cabeçalho formatado para a sequência
			codons = []  #Cria uma lista vazia para armazenar os códons da sequência atual

			#Iterando pela string da sequência, começando do "inicio_index" e pulando de 3 em 3 bases (quadro de leitura 1)
			for i in range(inicio_index, len(sequencia), 3):
				codon = sequencia[i:i+3]  #Extrai o slice de 3 bases
				if len(codon) == 3:  #Garante que o códon esteja completo (ignora "restos" no final da sequência)
					codons.append(codon)  #Adiciona o códon completo à lista

			#Juntando todos os códons da lista, separados por um espaço
			linha_codons_formatada = " ".join(codons)

			#Escrevendo a linha de códons formatada no arquivo
			arquivo_saida.write(f"{linha_codons_formatada}\n")

			#Adicionando uma linha em branco para separar as entradas no arquivo de saída (escolha estética)
			arquivo_saida.write("\n")

#Emitindo uma mensagem final para o usuário, indicando que o processo terminou
print("Processamento concluído. Verifique o arquivo \"Python_08.codons-3frames.nt\".")
