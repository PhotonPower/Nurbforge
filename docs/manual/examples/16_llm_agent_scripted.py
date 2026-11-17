"""The agent loop with a scripted stand-in for the model (no Ollama needed).

A real session replaces ScriptedClient by nurbforge.llm.OllamaClient; the agent only uses
its chat(messages, tools) method, which returns one assistant message per call.
"""
from nurbforge.llm import GeometryAgent


def call(tool, **args):
    return {"function": {"name": tool, "arguments": args}}


class ScriptedClient:
    """Replays prepared assistant messages and records what the agent sent."""

    def __init__(self, replies):
        self.replies = list(replies)
        self.requests = []

    def chat(self, messages, tools=None, temperature=0.0):
        self.requests.append(messages[-1])
        return self.replies.pop(0)


client = ScriptedClient([
    # 1st answer: the model creates an open shell by mistake (only 5 faces) ...
    {"role": "assistant", "content": "", "tool_calls": [
        call("make_rectangle", name="floor", origin=[0, 0, 0], edge_u=[1, 0, 0],
             edge_v=[0, 1, 0])]},
    {"role": "assistant", "content": "Created the floor plate."},
    # ... the verification report shows the result; the model repairs it with a solid
    {"role": "assistant", "content": "", "tool_calls": [
        call("solid", name="cube", kind="box", point=[0, 0, 0], size=[1, 1, 1]),
        call("check", name="cube", property="volume", expected=1)]},
    {"role": "assistant", "content": "Replaced the plate by a closed unit cube."},
    # second verification round: the model confirms without further tool calls
    {"role": "assistant", "content": "VERIFIED: closed unit cube 'cube', volume 1."},
])
agent = GeometryAgent(client, verify_rounds=2)
print("answer:", agent.ask("Make a closed unit cube at the origin."))
for name, args, result in agent.log:
    print(f"  [{name}] -> {result}")
print("last verification report sent to the model:")
print(agent.verification)
print("workspace:", sorted(agent.workspace.objects))
