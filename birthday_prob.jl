### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ b5743281-5981-48f8-bb2c-ac0845e71d8c
arr = [2,6,2,10]

# ╔═╡ 73c86261-52b8-4492-bd1e-f2045a5b9d9a
arr.==2
#Check if any element equal to 2

# ╔═╡ 79c6d135-b1df-49d0-90ce-60c33efd35fb
sum(arr.==2)
#retrun number of times 2 is there

# ╔═╡ 0868a708-99c8-410e-af62-f619e45a92fd
md"""
Using this method we can easilty write a code to check the probability.
"""

# ╔═╡ ee9dc083-2e00-4747-a09b-641b0665cd38
begin
	function prob(n::Int)
		favour = 0
		for ex_num=1:n
			months = rand(1:12,20)
			nums = [sum(months.==i) for i=1:12]
			if sum(nums.==3)==4 && sum(nums.==2)==4
				favour += 1
			end
		end
		probability = favour/n
	end
end

# ╔═╡ 6e077308-58d5-4a6b-9eae-88c06e0cce74
prob_10_000_000 = prob(10_000_000)

# ╔═╡ Cell order:
# ╠═b5743281-5981-48f8-bb2c-ac0845e71d8c
# ╠═73c86261-52b8-4492-bd1e-f2045a5b9d9a
# ╠═79c6d135-b1df-49d0-90ce-60c33efd35fb
# ╟─0868a708-99c8-410e-af62-f619e45a92fd
# ╠═ee9dc083-2e00-4747-a09b-641b0665cd38
# ╠═6e077308-58d5-4a6b-9eae-88c06e0cce74
