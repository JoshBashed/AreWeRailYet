import { RotateCcw, ZoomIn, ZoomOut } from 'lucide-react';
import { useEffect, useMemo, useRef, useState } from 'react';

type InteractiveMapProps = {
    src: string;
    alt: string;
};

const clamp = (value: number, min: number, max: number) =>
    Math.min(Math.max(value, min), max);

export default function InteractiveMap({ src, alt }: InteractiveMapProps) {
    const [scale, setScale] = useState(1);
    const [offset, setOffset] = useState({ x: 0, y: 0 });
    const surfaceRef = useRef<HTMLDivElement | null>(null);
    const dragState = useRef({
        active: false,
        originX: 0,
        originY: 0,
        startX: 0,
        startY: 0,
    });

    const transform = useMemo(
        () => `translate(${offset.x}px, ${offset.y}px) scale(${scale})`,
        [offset.x, offset.y, scale],
    );

    const setZoom = (next: number) => {
        setScale(clamp(next, 1, 3));
        if (next <= 1) {
            setOffset({ x: 0, y: 0 });
        }
    };

    const handlePointerDown = (event: React.PointerEvent<HTMLDivElement>) => {
        dragState.current.active = true;
        dragState.current.startX = event.clientX;
        dragState.current.startY = event.clientY;
        dragState.current.originX = offset.x;
        dragState.current.originY = offset.y;
        event.currentTarget.setPointerCapture(event.pointerId);
    };

    const handlePointerMove = (event: React.PointerEvent<HTMLDivElement>) => {
        if (!dragState.current.active) return;
        const dx = event.clientX - dragState.current.startX;
        const dy = event.clientY - dragState.current.startY;
        setOffset({
            x: dragState.current.originX + dx,
            y: dragState.current.originY + dy,
        });
    };

    const handlePointerUp = (event: React.PointerEvent<HTMLDivElement>) => {
        dragState.current.active = false;
        event.currentTarget.releasePointerCapture(event.pointerId);
    };

    const handleWheel = (event: React.WheelEvent<HTMLDivElement>) => {
        event.preventDefault();
        const zoomDelta = -event.deltaY / 300;
        setScale((current) => clamp(current + zoomDelta, 1, 3));
    };

    useEffect(() => {
        const surface = surfaceRef.current;
        if (!surface) return;

        const handleNativeWheel = (event: WheelEvent) => {
            event.preventDefault();
            const zoomDelta = -event.deltaY / 300;
            setScale((current) => clamp(current + zoomDelta, 1, 3));
        };

        surface.addEventListener('wheel', handleNativeWheel, {
            passive: false,
        });
        return () => surface.removeEventListener('wheel', handleNativeWheel);
    }, []);

    return (
        <div className='relative h-full w-full cursor-grab overflow-hidden rounded-3xl border border-zinc-800 bg-zinc-950 active:cursor-grabbing'>
            <div
                className='h-full w-full overflow-hidden'
                onPointerDown={handlePointerDown}
                onPointerLeave={handlePointerUp}
                onPointerMove={handlePointerMove}
                onPointerUp={handlePointerUp}
                onWheel={handleWheel}
                ref={surfaceRef}
                style={{ overscrollBehavior: 'contain', touchAction: 'none' }}
            >
                <img
                    alt={alt}
                    className='h-full w-full select-none object-contain'
                    draggable={false}
                    src={src}
                    style={{ transform, transformOrigin: 'center' }}
                />
            </div>
            <div className='absolute top-2 right-2 flex flex-col gap-1'>
                {[
                    {
                        action: () => setZoom(scale + 0.2),
                        icon: <ZoomIn size={16} />,
                        label: 'Zoom In',
                    },
                    {
                        action: () => setZoom(scale - 0.2),
                        icon: <ZoomOut size={16} />,
                        label: 'Zoom Out',
                    },
                    {
                        action: () => {
                            setZoom(1);
                            setOffset({ x: 0, y: 0 });
                        },
                        icon: <RotateCcw size={16} />,
                        label: 'Reset',
                    },
                ].map(({ label, icon, action }) => (
                    <button
                        aria-label={label}
                        className='flex h-8 w-8 cursor-pointer items-center justify-center rounded-full border border-zinc-800 bg-zinc-800/70 transition hover:bg-zinc-700/50'
                        key={label}
                        onClick={action}
                        title={label}
                        type='button'
                    >
                        {icon}
                    </button>
                ))}
            </div>
        </div>
    );
}
